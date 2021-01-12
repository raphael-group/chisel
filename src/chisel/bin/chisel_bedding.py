#!/usr/bin/env python2.7

import sys, os
import argparse

from multiprocessing import Lock, Value, Pool
from collections import Counter

import chisel

src = os.path.dirname(chisel.__file__)
from ..Utils import *
from chisel import Plotter


def parse_args():
    description = "CHISEL command to generate a BED file for each cell with the corresponding CHISEL's results."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='?', default='calls/calls.tsv', type=str, help="Input file with inferred copy numbers (default: calls/calls.tsv)")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("--rawcalls", required=False, default=False, action='store_true', help="Use raw copy numbers instead of consensus corrected ones (default: False)")
    parser.add_argument("--noextending", required=False, default=False, action='store_true', help="Merge consecutive bins only if they are neighboring (default: False, segments are extended to fill gaps)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError('ERROR: input file {} does not exist!'.format(args.INPUT))
    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")
    
    return {
        "input" : args.INPUT,
        "rundir" : args.rundir,
        "rawcalls" : args.rawcalls,
        "noextending" : args.noextending,
        "j" : args.jobs
    }


def main():
    log('Parsing and checking arguments', level='STEP')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Reading input', level='STEP')
    bins, pos, cells, iscorr = Plotter.read_cells(args['input'])
    log('Number of cells: {}'.format(len(cells)), level='INFO')
    log('Number of bins: {}'.format(len(pos)), level='INFO')

    log('Converting and writing BED files', level='STEP')
    make_beds(bins, pos, cells, iscorr, args)

    log('KTHXBYE', level='STEP')


def make_beds(bins, pos, cells, iscorr, args):
    jobs = (e for e in cells)
    bar = ProgressBar(total=len(cells), length=40, verbose=False)
    chk = (lambda g, e : bins[g][e]['CORR-CNS' if iscorr and not args['rawcalls'] else 'CNS'] if e in bins[g] else None)
    cns = {g : {e : chk(g, e) for e in bins[g]} for g in bins}
    initargs = (cns, pos, sum(g[2] - g[1] for g in pos), args['rundir'], args['noextending'])
    pool = Pool(processes=min(args['j'], len(cells)), initializer=init_making_bed, initargs=initargs)
    progress = (lambda e : bar.progress(advance=True, msg="Wrote cell {}".format(e)))
    res = map(progress, pool.imap_unordered(making_bed, jobs))
    pool.close()
    pool.join()


def init_making_bed(_cns, _pos, _totcov, _rundir, _noextending):
    global cns, pos, totcov, rundir, noextending
    cns = _cns
    pos = _pos
    totcov = _totcov
    rundir = _rundir
    noextending = _noextending


def making_bed(e):
    ecns = {g : cns[g][e] for g in pos}
    flag = {g : x == 0 or g[0] != pos[x-1][0] for x, g in enumerate(pos)}
    form = (lambda cns : '{}|{}'.format(*cns))
    
    with open(os.path.join(rundir, '{}.bed'.format(e)), 'w') as o:    
        def out(start, end):
            assert start[0] == end[0] and start[1] <= end[1] and start[2] <= end[2] and ecns[start] == ecns[end], 'Error on consecutive bins'
            if not noextending:
                o.write('\t'.join(map(str, [start[0], 0 if flag[start] else start[1], end[2], form(ecns[start])])) + '\n')
            else:
                o.write('\t'.join(map(str, [start[0], start[1], end[2], form(ecns[start])])) + '\n')
        
        start = pos[0]
        end = None
        precn = None
        tot = 0
        for x, g in enumerate(pos):
            if end is not None and (start[0] != g[0] or precn != ecns[g] or (noextending and end[2] != g[1])):
                out(start, end)
                tot += end[2] - start[1]
                start = g
                precn = ecns[g]
            end = g
        out(start, end)
        tot += end[2] - start[1]
        
    assert end == pos[-1], 'Error for the last bin'
    assert (not noextending or tot == totcov) and (noextending or tot >= totcov), 'Error in total length: {} written vs. {} expected'.format(tot, totcov)
    return e
    

if __name__ == '__main__':
    main()
