#!/usr/bin/env python2.7

import os, sys
import shlex
import argparse
import subprocess as sp

from multiprocessing import Lock, Value, Pool
from collections import Counter

from Utils import *


def parse_args():
    description = "Compute RDR from barcoded single-cell sequencing data."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-t","--tumor", required=True, type=str, help="Barcoded BAM file")
    parser.add_argument("-n","--normal", required=True, type=str, help="BAM file for matched normal sample")
    parser.add_argument("-b","--size", type=str, required=False, default="5Mb", help="Bin size, with or without \"kb\" or \"Mb\"")
    parser.add_argument("-s","--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-r","--reference", type=str, required=False, default="hg19", help="Name of the corresponding reference genome among \{hg18, hg19, hg38\} (default: hg19)")
    parser.add_argument("-m","--minreads", type=int, required=False, default=100000, help="Minimum number total reads to select cells (default: None)")
    parser.add_argument("-l","--cellslist", type=str, required=False, default=None, help="List of cells to select (default: None)")
    parser.add_argument("-c", "--chromosomes", type=str, required=False, default=' '.join(['chr{}'.format(i) for i in range(1, 23)]), help="Space-separeted list of chromosomes between apices (default: \"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\")")
    parser.add_argument("--outdir", required=False, default='./', type=str, help="Running directory where to write the list of selected cells (default: current directory)")
    args = parser.parse_args()

    if not os.path.isfile(args.tumor):
        raise ValueError("Specified tumor does not exist!")
    if not os.path.isfile(args.normal):
        raise ValueError("Specified normal does not exist!")
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome not found!")
    if not os.path.isfile(os.path.splitext(args.reference)[0] + '.dict'):
        raise ValueError("The dictionary .dict of the reference genome not found! Please remember to index it.")
    if not os.path.isdir(args.outdir):
        raise ValueError("Running directory does not exists: {}".format(args.outdir))

    size = 0
    try:
        if args.size[-2:] == "kb":
            size = int(args.size[:-2]) * 1000
        elif args.size[-2:] == "Mb":
            size = int(args.size[:-2]) * 1000000
        else:
            size = int(args.size)
    except:
	raise ValueError("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!")

    samtools = args.samtools
    if not samtools:
        samtools = "samtools"
    if which(samtools) is None:
        raise ValueError("samtools has not been found or is not executable!")
    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")
    if args.minreads is not None and args.cellslist is not None:
        raise ValueError("Only one between number of reads or list of cells can be specified!")
    if args.minreads and args.minreads < 1:
        raise ValueError("Minimum number of reads must be positive!")
    if args.cellslist and not os.path.isfile(args.cellslist):
        raise ValueError("Cell list does not exist!")

    return {
        'tumor' : args.tumor,
        'normal' : args.normal,
        'bins' : size,
        'samtools' : samtools,
        'J' : args.jobs,
        'ref' : args.reference,
        'minreads' : args.minreads,
        'list' : args.cellslist,
        'chrs' : args.chromosomes.split(),
        'outdir' : args.outdir
    }


def main():
    log('Parsing and checking arguments')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['{} : {}'.format(a, args[a]) for a in args]), level='INFO')

    log('Computing bins')
    bins = get_bins(ref=args['ref'], chrs=args['chrs'], bsize=args['bins'])

    log('Counting reads on normal')
    counts = counting_normal(args['normal'], bins, args['samtools'], args['J'])

    log('Counting reads on barcoded cells')
    counts = counting_cells(counts, args['tumor'], bins, args['samtools'], args['J'])

    log('Evaluating set of found cells')
    if args['list'] is None:
        names = set(e for c in counts for b in counts[c] for e in counts[c][b])
        cells = names - {'normal'}
    else:
        clist = set()
        with open(args['list'], 'r') as i:
            for l in i:
                clist.add(l.strip().replace(',','\t').split()[0].replace('-1', ''))
        names = set(e for c in counts for b in counts[c] for e in counts[c][b] if e in clist)

    log('Computing total numbers of sequenced reads')
    total = reduce(inupdate, (Counter(counts[c][b]) for c in counts for b in counts[c]))

    log('Selecting cells')
    if args['minreads']:
        names = set(e for e in total if total[e] >= args['minreads'])
        cells = names - {'normal'}
    log('Number of selected cells: {}'.format(len(cells)), level='INFO')

    ftot = os.path.join(args['outdir'], 'total.tsv')
    log('Writing the totals in {}'.format(ftot), level='INFO')
    with open(ftot, 'w') as o:
        o.write('{}\t{}\n'.format('normal', total['normal']))
        o.write('\n'.join(['{}\t{}'.format(e, total[e]) for e in cells]))

    log('Estimating RDR')
    scale = {e : float(total['normal']) / float(total[e]) for e in cells}
    ratio = (lambda c, b, e : (float(counts[c][b][e]) / float(counts[c][b]['normal'])) if counts[c][b]['normal'] > 0 else 0.0)
    rec = (lambda c, b, e, rdr : '{}\t{}\{}\t{}'.format(c, b, e, rdr))

    for c in sorted(counts, key=orderchrs):
        for b in sorted(counts[c], key=(lambda x : x[0])):
            for e in sorted(set(counts[c][b].keys()) & cells):
                print '\t'.join(map(str, [c, b[0], b[1], e, counts[c][b]['normal'], counts[c][b][e], ratio(c, b, e) * scale[e]]))

    log('KTHXBYE')


def counting_normal(normal, bins, samtools, J):
    jobs = [(c, b) for c in bins for b in bins[c]]
    lock = Lock()
    counter = Value('i', 0)

    initargs = (lock, counter, len(jobs), normal, samtools)
    pool = Pool(processes=min(J, len(jobs)), initializer=init_extracting_normal, initargs=initargs)

    counts = defaultdict(lambda : defaultdict(lambda : dict()))
    try:
        for c, b, rd in pool.imap_unordered(extracting_normal, jobs):
            assert 'normal' not in counts[c][b]
            if rd != '':
                counts[c][b]['normal'] = int(rd.strip())
            else:
                counts[c][b]['normal'] = 0
    except Exception, e:
        pool.close()
        pool.terminate()
        raise RuntimeError("ERROR: " + str(e))
        sys.exit(1)

    return counts


def init_extracting_normal(lock, counter, _l, _normal, sam):
    global bar, cmd_sam
    bar = ProgressBar(total=_l, length=40, lock=lock, counter=counter, verbose=False)
    cmd_sam = "{} view -F 1796 -q 13 -c {} {}:{}-{}".format(sam, _normal, "{}", "{}", "{}")


def extracting_normal(job):
    c, b = job
    cmd = cmd_sam.format(c, b[0], b[1])
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    bar.progress(advance=True, msg="Counting normal on {}:{}-{}".format(c, b[0], b[1]))
    return (c, b, stdout.strip())


def counting_cells(counts, tumor, bins, samtools, J):
    jobs = [(c, b) for c in bins for b in bins[c]]
    bar = ProgressBar(total=len(jobs), length=40, verbose=False)

    initargs = (tumor, samtools)
    pool = Pool(processes=min(J, len(jobs)), initializer=init_extracting, initargs=initargs)

    for c, b, rd in pool.imap_unordered(extracting, jobs):
        if rd != '':
            for l in rd.strip().split('\n'):
                p = l.split()
                assert p[0] not in counts[c][b]
                counts[c][b][p[0]] = int(p[1])
        bar.progress(advance=True, msg="Extracted barcodes on {}:{}-{}".format(c, b[0], b[1]))

    return counts


def init_extracting(_tumor, sam):
    global cmd_sam, cmd_awk
    cmd_sam = "{} view -F 1796 -q 13 {} {}:{}-{}".format(sam, _tumor, "{}", "{}", "{}")
    cmd_awk = shlex.split("awk 'BEGIN{} { if(match($0, /CB:Z:[ACGT]+/)) { X[substr($0, RSTART+5, RLENGTH-5)]++ } } END{ for(i in X) print i, X[i] }'")


def extracting(job):
    c, b = job
    cmd = cmd_sam.format(c, b[0], b[1])
    sam = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(cmd_awk, stdin=sam.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return (c, b, stdout.strip())


def get_bins(ref, chrs, bsize):
    ends = {}
    ends["hg38"] = [-1, 248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]
    ref = "hg38"

    ln = (lambda c : ends[ref][int(''.join([l for l in c if l.isdigit()]))])
    fl = (lambda l, c : l + [ln(c)] if l[-1] < ln(c) else (l if l[-1] == ln(c) else []))
    bk = (lambda c : fl(list(range(0, ln(c), bsize)), c))
    bins = {c : sorted(zip(bk(c)[:-1], bk(c)[1:]), key=(lambda x : x[0])) for c in chrs}
    assert False not in set(len(bins[c]) > 0 for c in bins), "Binning failed for some chromosomes"

    return bins


# def get_bins(ref, chrs, bsize):
#     ends = {c : None for c in chrs}
#     refdict = os.path.splitext(ref)[0] + '.dict'
#     with open(refdict, 'r') as i:
#         for l in i:
#             if '@SQ' in l:
#                 assert 'SN:' in l and 'LN:' in l
#                 c = l.split('SN:')[1].split()[0]
#                 if c in chrs:
#                     end = int(l.split('LN:')[1].split()[0])
#                     ends[c] = end
#     if None in ends.values():
#         err("The following chromosomes have not been found in the dictionary of the reference genome: \n\t{}".format(','.join([c for c in ends if ends[c] == None])))
#     print ends

#     fl = (lambda l, c : l + [ends[c]] if l[-1] < ends[c] else (l if l[-1] == ends[c] else []))
#     bk = (lambda c : fl(list(range(0, ends[c], bsize)), c))
#     bins = {c : sorted(zip(bk(c)[:-1], bk(c)[1:]), key=(lambda x : x[0])) for c in chrs}
#     assert False not in set(len(bins[c]) > 0 for c in bins), "Binning failed for some chromosomes"

#     return bins


if __name__ == '__main__':
    main()
