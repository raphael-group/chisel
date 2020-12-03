#!/usr/bin/env python2.7

import os, sys
import argparse
import bisect
import math
import multiprocessing as mp

from multiprocessing import Lock, Value, Pool
from collections import defaultdict
from collections import Counter
from collections import deque
from functools import reduce

import numpy as np
import scipy.stats

from Utils import *


def parse_args(args):
    description = "Compute RDR from barcoded single-cell sequencing data."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-r","--rdr", required=True, type=str, help="RDR file")
    parser.add_argument("-b","--baf", required=True, type=str, help="BAF file")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-k", "--blocksize", required=False, type=str, default="50kb", help="Size of the haplotype blocks (default: 50kb, use 0 to disable)")
    parser.add_argument("-a", "--significance", required=False, type=float, default=0.05, help="Significance level to estimate the maximum shift error for BAF = 0.5 (default: 0.05)")
    parser.add_argument("-q", "--restarts", required=False, type=int, default=100, help="Number of restarts for EM (default: 100)")
    parser.add_argument("-t", "--bootstrap", required=False, type=int, default=100, help="Number of bootastrapping points to estimate maximum shift error (default: 100)")
    parser.add_argument("-e", "--maxerror", required=False, type=float, default=None, help="Maximum shift error for identification of BAF = 0.5 (default: not used, when specified the value is used instead of estimation)")
    parser.add_argument("-E", "--minerror", required=False, type=float, default=0.001, help="Minimum shift error for identification of BAF = 0.5 (default: 0.001)")
    parser.add_argument("-l", "--listofcells", required=False, type=str, default=None, help="List of cells to include (default: None)")
    parser.add_argument("-s", "--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    args = parser.parse_args(args)

    if not os.path.isfile(args.rdr):
        raise ValueError("RDR file does not exist!")
    if not os.path.isfile(args.baf):
        raise ValueError("BAF file does not exist!")

    if args.jobs == 0:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")
    if args.restarts < 1:
        raise ValueError("The number of restarts must be positive!")
    if args.bootstrap < 1:
        raise ValueError("The number of bootstrapping points must be positive!")
    if not 0.0 <= args.significance <= 1.0:
        raise ValueError("The maxerror must be in [0, 1]!")
    if args.maxerror is not None and not 0.0 <= args.maxerror <= 0.5:
        raise ValueError("The maxerror must be in [0, 0.5]!")
    if args.minerror is not None and not 0.0 <= args.minerror <= 0.5:
        raise ValueError("The minerror must be in [0, 0.5]!")
    if args.listofcells is not None and not os.path.isfile(args.listofcells):
        raise ValueError("The list of cells does not exist!")
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")

    blocksize = 0
    try:
        if args.blocksize[-2:] == "kb":
            blocksize = int(args.blocksize[:-2]) * 1000
        elif args.blocksize[-2:] == "Mb":
            blocksize = int(args.blocksize[:-2]) * 1000000
        else:
            blocksize = int(args.blocksize)
    except:
        raise ValueError("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!")

    if blocksize == 0:
        blocksize = None   

    return {
        'rdr' : args.rdr,
        'baf' : args.baf,
        'j' : args.jobs,
        'blocksize' : blocksize,
        'restarts' : args.restarts,
        'bootstrap' : args.bootstrap,
        'significance' : args.significance,
        'maxerror' : args.maxerror,
        'minerror' : args.minerror,
        'listofcells' : args.listofcells,
        'seed' : args.seed
    }


def main(args=None, stdout_file=None):
    log('Parsing and checking arguments')
    args = parse_args(args)
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]), level='INFO')

    cells = None
    if args['listofcells'] is not None:
        log('Read list of cells')
        cells = read_listofcells(args['listofcells'])

    log('Reading RDR')
    rdr = read_rdr(args['rdr'], cells=cells)

    log('Reading BAF')
    cA, cB, bulk = read_baf(args['baf'], cells=cells)

    log('Combining')
    rb = combo(rdr, cA, cB, bulk, args)

    if stdout_file is not None:
        stdout_f = open(stdout_file, 'w')

    log('Printing combined RDR and BAF')
    baf = (lambda A, B : (float(B) / float(A+B)) if A+B > 0 else 0.5)
    gen = ((c, b, e) for c in sorted(rb, key=orderchrs) for b in sorted(rb[c], key=(lambda x : x[0])) for e in sorted(rb[c][b]))
    for c, b, e in gen:
        formb = (lambda c, b, e : [rb[c][b][e]['BAF'][0], rb[c][b][e]['BAF'][1], baf(*rb[c][b][e]['BAF'])])
        formr = (lambda c, b, e : [rb[c][b][e]['normalcount'], rb[c][b][e]['readcount'], rb[c][b][e]['RDR']])
        line = '\t'.join(map(str, [c, b[0], b[1], e] + formr(c, b, e) + formb(c, b, e)))
        if stdout_file is not None:
            stdout_f.write(line + '\n')
        else:
            print line

    if stdout_file is not None:
        stdout_f.close()


def read_rdr(f, cells):
    d = (lambda : {'RDR' : 0.0, 'readcount' : 0, 'normalcount' : 0})
    rdr = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : d())))
    with open(f, 'r') as i:
        for l in (l for l in i if l[0] != '#' and len(l) > 1):
            p = l.strip().split()
            if cells is None or p[3] in cells:
                assert len(p) == 7 and p[3] not in rdr[p[0]][(int(p[1]), int(p[2]))]
                rdr[p[0]][(int(p[1]), int(p[2]))][p[3]]['normalcount'] = int(p[-3])
                rdr[p[0]][(int(p[1]), int(p[2]))][p[3]]['readcount'] = int(p[-2])
                rdr[p[0]][(int(p[1]), int(p[2]))][p[3]]['RDR'] = float(p[-1])
    return {c : dict(rdr[c]) for c in rdr}


def read_baf(f, cells):
    countA = defaultdict(lambda : defaultdict(lambda : dict()))
    countB = defaultdict(lambda : defaultdict(lambda : dict()))
    bulk = defaultdict(lambda : defaultdict(lambda : (0, 0)))
    form = (lambda p0, p1, p2, p3, p4 : (p0, int(p1), p2, int(p3), int(p4)))
    with open(f, 'r') as i:
        for l in (l for l in i if l[0] != '#' and len(l) > 1):
            c, o, e, a, b = form(*l.strip().split())
            if (a + b) > 0 and (cells is None or e in cells):
                assert e not in countA[c][o] and e not in countB[c][o]
                countA[c][o][e] = a
                countB[c][o][e] = b
                bulk[c][o] = (bulk[c][o][0] + a, bulk[c][o][1] + b)
    return {c : dict(countA[c]) for c in countA}, {c : dict(countB[c]) for c in countB}, {c : dict(bulk[c]) for c in bulk}


def read_listofcells(f):
    cells = set()
    with open(f, 'r') as i:
        for l in (l for l in i if l[0] != '#' and len(l) > 1):
            cells.add(l.strip().split()[0])
    return cells
            

def combo(rdr, cA, cB, bulk, args):
    np.random.seed(args['seed'])
    rb = defaultdict(lambda : dict())
    cells = set(e for c in rdr for b in rdr[c] for e in rdr[c][b])
    jobs = ((c, b, np.random.randint(1000)) for c in rdr for b in rdr[c])
    njobs = sum(len(rdr[c].keys()) for c in rdr)
    snps = {c : sorted(cA[c].keys()) if c in cA else [] for c in rdr}
    bar = ProgressBar(total=njobs, length=40, verbose=False)

    initargs = (snps, cA, cB, bulk, args)
    pool = Pool(processes=min(args['j'], njobs), initializer=init, initargs=initargs)
    counts = (lambda c, b, e : rdr[c][b][e].items())
    for c, b, A, B in pool.imap_unordered(combine, jobs):
        rb[c][b] = {e : dict(counts(c, b, e) + [('BAF', (A[e], B[e]) if e in A else (0, 0))]) for e in cells}
        bar.progress(advance=True, msg="Combined bin {}:{}-{}".format(c, b[0], b[1]))

    return rb
    

def init(_snps, _cA, _cB, _bulk, args):
    global snps, cA, cB, bulk, blocksize, restarts, boot, alpha, maxerror, minerror
    snps = _snps
    cA = _cA
    cB = _cB
    bulk = _bulk
    blocksize = args['blocksize']
    restarts = args['restarts']
    boot = args['bootstrap']
    alpha = args['significance']
    maxerror = args['maxerror']
    minerror = args['minerror']


def combine(job):
    c, b, seed = job
    np.random.seed(seed)
    L = bisect.bisect_left(snps[c], b[0])
    R = bisect.bisect_right(snps[c], b[1])

    if L >= R:
        return c, b, dict(), dict()

    snpsLR = snps[c][L:R]
    assert all(snpsLR[x - 1] < snpsLR[x] if x > 0 else True for x, o in enumerate(snpsLR))

    if blocksize:
        que = deque(snpsLR)
        assert sorted(snpsLR) == list(que) and b[0] <= que[0] and que[-1] <= b[1]
        omap = {}
        blocks = {}
        for bk in range(b[0], b[1]+1, blocksize):
            block = (0, 0)
            while que and bk <= que[0] < bk + blocksize:
                o = que.popleft()
                block = (block[0] + bulk[c][o][0], block[1] + bulk[c][o][1])
                omap[o] = bk
            if sum(block) > 0:
                blocks[bk] = block
        assert set(omap.values()) == set(blocks.keys())
        assert set(snpsLR) == set(omap.keys())
                
        allblocks = blocks.values()
        nis = [sum(block) for block in allblocks]
        xis = [block[0] if np.random.random() < 0.5 else block[1] for block in allblocks]
        beta = max((EM(ns=nis, xs=xis, start=np.random.randint(low=1, high=49)/100.0) for r in xrange(restarts)), key=(lambda x : x[1]))[0]
        if maxerror is None:
            thres = max(minerror, est_error(ns=nis, significance=alpha, restarts=restarts, bootstrap=boot))
        else:
            thres = maxerror
            
        if thres <= abs(beta - 0.5): # <= 0.25:
            minel = (lambda p : 0 if p[0] < p[1] else 1)
            sumpairs = (lambda I : reduce((lambda x, y : (x[0] + y[0], x[1] + y[1])), I))
            if minel(sumpairs(allblocks)) == 0:
                swap = {o : False if blocks[omap[o]][0] < blocks[omap[o]][1] else True for o in snpsLR}
            else:
                swap = {o : False if blocks[omap[o]][1] < blocks[omap[o]][0] else True for o in snpsLR}
        else:
            bkswap = {bk : False if np.random.random() < 0.5 else True for bk in blocks}
            swap = {o : True if bkswap[omap[o]] else False for o in snpsLR}
    else:
        swap = {o : False for o in snpsLR}
                
    A = reduce(inupdate, (Counter(cA[c][o] if not swap[o] else cB[c][o]) for o in snpsLR))
    B = reduce(inupdate, (Counter(cB[c][o] if not swap[o] else cA[c][o]) for o in snpsLR))
    return c, b, A, B


def EM(ns, xs, start, tol=10**-6):
    nis = np.array(ns)
    xis = np.array(xs)
    assert nis.size == xis.size and 0 < start < 1 and np.all(nis >= xis)
    if np.all(np.logical_or(nis == xis, xis == 0)):
        return 0.0, 0.0
    else:
        beta = start
        prev = None
        while prev is None or abs(prev - beta) >= tol:
            prev = beta
            #print 'E-step'
            assert 0 + tol < beta < 1 - tol, (beta, nis, xis, start)
            M = (nis - 2*xis) * np.log(beta) + (2*xis - nis) * np.log(1.0 - beta)
            M = np.exp(np.clip(a=M, a_min=-100, a_max=100))
            his = np.reciprocal(1 + M)
            #print 'M-step'
            beta = float(np.sum(nis * (1 - his) + xis * (2*his - 1))) / float(np.sum(nis))
            #print 'BETA = {}'.format(beta)
        assert 0 + tol < beta < 1 - tol, (beta, nis, xis, start)
        lpmf = scipy.stats.binom.logpmf
        loglh = float(np.sum(his * lpmf(k=xis, n=nis, p=beta) + (1 - his) * lpmf(k=xis, n=nis, p=1-beta)))
        return beta, loglh


def est_error(ns, significance=0.05, restarts=50, bootstrap=100):
    nis = np.array(ns)
    N = nis.size
    genstart = (lambda : np.random.randint(low=1, high=49)/100.0)
    runEM = (lambda xis : max((EM(ns=nis, xs=xis, start=genstart()) for x in xrange(restarts)), key=(lambda x : x[1]))[0])
    mirror = (lambda b : min(b, 1 - b))
    genneu = scipy.stats.binom.rvs
    betas = sorted(mirror(runEM(genneu(n=nis, p=0.5, size=len(nis)))) for x in xrange(bootstrap))
    betas = betas[int(round(len(betas) * significance)):]
    return 0.5 - betas[0]

    
if __name__ == '__main__':
    main()

