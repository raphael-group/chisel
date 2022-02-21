#!/usr/bin/env python2.7

import os
import argparse
import shlex, shutil
import multiprocessing as mp

from multiprocessing import Lock, Value, Pool
from collections import defaultdict
from collections import Counter

import numpy as np

import chisel

src = os.path.dirname(chisel.__file__)
from ..Utils import *
from ..RDREstimator import *


def parse_args():
    description = "CHISEL command to generate a pseudo-matched normal sample by extracting diploid cells from a barcoded single-cell BAM file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", type=str, help="Barcoded single-cell BAM file")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-e","--threshold", type=float, required=False, default=0.9, help="Minimum fraction of diploid genome to select diploid cells (default: 0.9)")
    parser.add_argument("-b","--size", type=str, required=False, default="5Mb", help="Bin size, with or without \"kb\" or \"Mb\"")
    parser.add_argument("-c", "--chromosomes", type=str, required=False, default=' '.join(['chr{}'.format(i) for i in range(1, 23)]), help="Space-separeted list of chromosomes between apices (default: \"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\")")
    parser.add_argument("-m","--minreads", type=int, required=False, default=100000, help="Minimum number total reads to select cells (default: 100000)")
    parser.add_argument("--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("--tmpdir", required=False, default='_TMP', type=str, help="Temporary directory in running directory (default: _TMP)")
    parser.add_argument("-n","--normal", required=False, type=str, default="pseudonormal.bam", help="Name of the generated pseudo matched-normal BAM file (default: pseudonormal.bam)")
    args = parser.parse_args()

    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    tmpdir = os.path.join(args.rundir, args.tmpdir)
    if os.path.isdir(tmpdir):
        raise ValueError("Temporary directory already exists within specified running direcotyr: {}".format(tmpdir))    
    if not os.path.isfile(args.INPUT):
        raise ValueError("Barcoded single-cell BAM file does not exist: {}".format(args.INPUT))
    if not args.normal[:-4] != ".bam":
        raise ValueError("The provided output name does not end in .bam: {}".format(args.normal))
    if not (0.0 <= args.threshold <= 1.0):
        raise ValueError("The provided threshold is not in [0, 1]: {}".format(args.threshold))
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome file does not exist: {}".format(args.reference))
    if args.minreads < 1:
        raise ValueError("The minimum number of reads must be positive!")

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

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")

    samtools = args.samtools
    if not samtools:
        samtools = "samtools"
    if which(samtools) is None:
        raise ValueError("samtools has not been found or is not executable!")

    return {
        "rundir" : args.rundir,
        "tmpdir" : tmpdir,
        "tumor" : args.INPUT,
        "thres" : args.threshold,
        "normal" : args.normal,
        "reference" : args.reference,
        "binsize" : size,
        "chromosomes" : args.chromosomes,
        "minreads" : args.minreads,
        "samtools" : samtools,
        "jobs" : args.jobs
    }


def main():
    log('Parsing and checking arguments')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Computing bins')
    chrs = args['chromosomes'].split()
    bins = get_bins(args['reference'], chrs, args['binsize'], bams=[args['tumor']], samtools=args['samtools'])
    counts = defaultdict(lambda : defaultdict(lambda : dict()))
    assert len(chrs) == len(bins.keys())
    chrs = bins.keys()
    
    log('Counting reads on barcoded cells')
    counts = counting_cells(counts, args['tumor'], bins, args['samtools'], args['jobs'], args['prefix'], args['suffix'])
    cells = set(e for c in counts for b in counts[c] for e in counts[c][b])
    
    log('Computing total numbers of sequenced reads')
    total = reduce(inupdate, (Counter(counts[c][b]) for c in counts for b in counts[c]))

    log('Selecting all cells to consider for the analysis')
    if args['minreads']:
        cells = set(e for e in total if total[e] >= args['minreads'])
    log('Number of selected cells: {}'.format(len(cells)), level='INFO')

    log('Selecting diploid cells')
    diploid = sorted(set(filter((lambda e : isdiploid(counts, e, args['thres'])), cells)))
    dlist = os.path.join(args['rundir'], 'diploid.tsv')
    with open(dlist, 'w') as o:
        o.write('\n'.join(diploid) + '\n')
    log('Number of identified diploid cells: {}'.format(len(diploid)), level='INFO')
    cov = (float(sum(total[e] for e in diploid)) * 100.0) / float(sum(b[1] - b[0] for c in counts for b in counts[c]))
    log('Approximate sequencing coverage of pseudo matched-normal sample: {}'.format(cov), level='INFO')

    log('Extracting sequencing reads from selected diploid cells')
    extracting_diploid(args['tumor'], args['samtools'], chrs, args['tmpdir'], dlist, args['jobs'])

    log('Merging and extracted sequencing reads and indexing the output pseduo matched-normal sample')
    merging_diploid(chrs, args['tmpdir'], args['samtools'], os.path.join(args['rundir'], args['normal']))

    log('Removing temporary files')
    shutil.rmtree(args['tmpdir'])

    log('KTHXBYE')


def isdiploid(counts, cell, THRES):
    rdr = np.array([counts[c][b][cell] for c in counts for b in counts[c] if cell in counts[c][b] and counts[c][b][cell] > 0])
    base = np.sum(rdr) / float(rdr.shape[0])
    assert base > 0, "Found a cell with no sequencing reads"
    rdr = rdr / base
    avg = 2.0 / (np.sum(rdr) / float(rdr.shape[0]))
    dip = (lambda t : np.sum(np.rint(t * rdr) == 2))
    scale = max((avg + (x/100.0)*d for x in xrange(0, 100+1, 1) for d in {-1, 1}), key=dip)
    return (dip(scale) / float(rdr.shape[0])) >= THRES


def extracting_diploid(bam, samt, chrs, tmpdir, dlist, J):
    lock = Lock()
    counter = Value('i', 0)
    assert not os.path.isdir(tmpdir)
    os.mkdir(tmpdir)
    initargs = (lock, counter, len(chrs), bam, samt, tmpdir, dlist)
    pool = Pool(processes=min(J, len(chrs)), initializer=init_extracting_diploid, initargs=initargs)
    res = {o for o in pool.imap_unordered(extract, chrs)}
        #if o.strip() != '':
            #raise ValueError("SAMtools raised the following error during extraction ofsequencing reads: {}".format(o))
    return

        
def init_extracting_diploid(lock, counter, l, bam, samt, _tmpdir, dlist):
    global bar, cmd_sam, cmd_gre, cmd_com, tmpdir
    bar = ProgressBar(total=l, length=min(l, 40), lock=lock, counter=counter, verbose=False)
    cmd_sam = "{} view -h -F 1796 -q 13 {} {}".format(samt, bam, "{}")
    cmd_gre = "grep -F -f {} -e \"@HD\" -e \"@SQ\" -e \"@RG\" -e \"@PG\" -e \"@CO\"".format(dlist)
    cmd_com = "{} sort -O bam -o {} -T {}".format(samt, "{}", "{}")
    tmpdir = _tmpdir


def extract(c):
    cmd = cmd_sam.format(c)
    out = os.path.join(tmpdir, '{}.bam'.format(c))
    tmp = os.path.join(tmpdir, '_TMP_{}'.format(c))
    os.mkdir(tmp)
    sam = sp.Popen(shlex.split(cmd_sam.format(c)), stdout=sp.PIPE, stderr=sp.PIPE)
    gre = sp.Popen(shlex.split(cmd_gre), stdin=sam.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(shlex.split(cmd_com.format(out, tmp)), stdin=gre.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return stderr.strip()


def merging_diploid(chrs, tmpdir, samt, out):
    cfiles = map((lambda c : os.path.join(tmpdir, '{}.bam'.format(c))), sorted(chrs, key=orderchrs))
    assert all(os.path.isfile(f) for f in cfiles), "Extracted reads are missing for some files!"
    cmd = "{} merge -f {} {}".format(samt, out, ' '.join(cfiles))
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    #if stderr.strip() != '':
    #    raise ValueError("SAMtools merging terminated with the following error: {}".format(stderr))
    cmd = "{} index {}".format(samt, out)
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    #if stderr.strip() != '':
    #    raise ValueError("SAMtools indexing terminated with the following error: {}".format(stderr))
    return


if __name__ == '__main__':
    main()
