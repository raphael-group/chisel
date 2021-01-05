#!/usr/bin/env python2.7

import os, sys
import shlex
import argparse
import subprocess as sp

from multiprocessing import Lock, Value, Pool
from collections import defaultdict

from Utils import *


def parse_args():
    description = "Cell-specific allele counting for a given list of point mutations." #which must be provided as a stdin stream (symbol '-' must be used in this case) or as the name of a file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-l","--listmutations", type=str, required=True, help="List of TSV phased genomic positions (TAB-seprated format '#CHR POS REF VAR')")
    parser.add_argument("-t","--tumor", required=True, type=str, help="BAM file for matched normal sample")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-s","--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-c","--listcells", type=str, required=False, default=None, help="File where first column contains all the cells to consider (default: not used)")
    args = parser.parse_args()

    if not os.path.isfile(args.tumor):
        raise ValueError("Specified tumor does not exist!")
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome does not exist!")
    if args.listcells is not None and not os.path.isfile(args.listcells):
        raise ValueError("Specified list of cells does not exist!")

    samtools = args.samtools
    if not samtools:
        samtools = "samtools"
    if which(samtools) is None:
        raise ValueError("samtools has not been found or is not executable!")

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")

    return {
        'tumor' : args.tumor,
        'mutations' : args.listmutations,
        'ref' : args.reference,
        'samtools' : args.samtools,
        'J' : args.jobs,
        'list' : args.listcells
    }


def main():
    log('Parsing and checking arguments')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]), level='INFO')

    log('Extracting genomic positions of given mutations')
    mutations = read_mutations(args['mutations'])
    log('Chromosomes analyzed: {}'.format(','.join(sorted(mutations, key=orderchrs))), level='INFO')
    log('Total number of given mutations: {}'.format(sum(len(mutations[c]) for c in mutations)), level='INFO')

    log('Extracting allele counts of mutations for all cells')
    amut = extracting(args, mutations)

    if args['list']:
        log('Reading cell list')
        with open(args['list'], 'r') as i:
            cells = set(l.strip().split()[0].replace('-1', '') for l in i if len(l) > 1 and l[0] != '#')

    log('Writing A/B counts for selected phased SNPs across selected cells')
    print '\t'.join(['#CHR', 'POS', 'CELL', 'MUT', 'MUTCOV', 'COV'])
    for c, o, e in ((c, o, e) for c in sorted(amut, key=orderchrs) for o in sorted(amut[c]) for e in sorted(amut[c][o])):
        print '\t'.join(map(str, [c, o, e, mutations[c][o], amut[c][o][e][mutations[c][o]], sum(amut[c][o][e].values())]))

    log('KTHXBYE')


def read_mutations(f):
    mutations = defaultdict(lambda : dict())
    chrs = map(str, range(1, 23))
    with open(f, 'r') as i:
        for l in i:
            if len(l) > 1 and l[0] != '#':
                p = l.strip().split()
                c = p[0]
                if ''.join([l for l in c if l.isdigit()]) not in chrs:
                    continue
                try:
                    o = int(p[1])
                    v = p[3]
                    assert o not in mutations[c]
                    if v[0] in {'A', 'C', 'G', 'T'} or v[0] in {'+', '-'}:
                        mutations[c][o] = 'N' if v[0] == '-' else (v[1] if v[0] == '+' else v[0])
                        assert mutations[c][o] in {'A', 'C', 'G', 'T', 'N'} 
                except ValueError:
                    pass
    return mutations


def extracting(args, mutations):
    jobs = ((c, o) for c in mutations for o in mutations[c])
    njobs = sum(len(mutations[c]) for c in mutations)
    countawk = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'count.awk')
    bar = ProgressBar(total=njobs, length=40, verbose=False)
    
    initargs = (args['tumor'], args['samtools'], countawk)
    pool = Pool(processes=min(args['J'], njobs), initializer=init_extracting, initargs=initargs)

    ACGT = (lambda : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'N' : 0})
    amut = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : ACGT())))
    amut = {c : {o : defaultdict(lambda : ACGT()) for o in mutations[c]} for c in mutations}
    for c, o, l in pool.imap_unordered(counting_cell, jobs):
        if l != '':
            for a in l.strip().split('\n'):
                e, al, count = tuple(a.split())
                amut[c][o][e][al] += int(count)
        bar.progress(advance=True, msg="Extracted SNP {}:{}".format(c, o))   

    return {c : {o : dict(filter(lambda (e, al) : sum(al.values()) > 0, amut[c][o].items())) for o in amut[c]} for c in amut}
    

def init_extracting(_tumor, _sam, countawk):
    global cmd_sam, cmd_awk
    cmd_sam = "{} view -F 1796 -q 13 {} {}:{}-{}".format(_sam, _tumor, '{}', '{}', '{}')
    cmd_awk = 'awk -v TAG="{}" -f {}'.format('{}', countawk)


def counting_cell(job):
    sam = sp.Popen(shlex.split(cmd_sam.format(job[0], job[1], job[1])), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(shlex.split(cmd_awk.format(job[1])), stdin=sam.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return (job[0], job[1], stdout)


if __name__ == '__main__':
    main()
