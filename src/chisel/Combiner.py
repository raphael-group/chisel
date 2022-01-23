#!/usr/bin/env python2.7

import os, sys
import argparse
import bisect
import math
import multiprocessing as mp
import warnings

from multiprocessing import Lock, Value, Pool
from collections import defaultdict
from collections import Counter
from collections import deque
from functools import reduce

import numpy as np
import scipy.stats
import statsmodels.api as sm

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
    parser.add_argument("--minimumsnps", required=False, type=float, default=0.0, help="Minimum SNP density per kb (default: 0.0, reasonable values are 0.01, 0.02, etc.)")
    parser.add_argument("--missingsnps", required=False, type=str, default=None, help="A,B counts for genomic bins without minimum minimum SNP density (default: 0,0 i.e. BAF=0.5)")
    parser.add_argument("--gccorr", required=False, type=str, default=None, help="The reference genome to apply additional GC correction in addition to using a matched-normal sample (default: None)")
    parser.add_argument("--nophasecorr", required=False, default=False, action='store_true', help="Disable correction for given phasing bias (default: enabled)")
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
        raise ValueError("The random seed must be positive!")
    if args.gccorr is not None and not os.path.isfile(args.gccorr):
        raise ValueError("The specified reference genome for additional GC correction does not exist!")
    if args.minimumsnps < 0.0:
        raise ValueError("The minimum SNP density must be >= 0.0!")

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

    if args.missingsnps is not None:
        if ',' not in args.missingsnps:
            raise ValueError("missingsnps parameter has the wrong format! Please specify A,B with A and B integers.")
        val = args.missingsnps.split(',')
        if len(val) != 2 and any(not l.isdigit() for l in val[0]) or any(not l.isdigit() for l in val[1]):
            raise ValueError("missingsnps parameter has the wrong format! Please specify A,B with A and B integers.")
        missingsnps = (int(val[0]), int(val[1]))
    else:
        missingsnps = None

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
        'seed' : args.seed,
        'minimumsnps' : args.minimumsnps,
        'missingsnps' : missingsnps,
        'gccorr' : args.gccorr,
        "phasecorr" : not args.nophasecorr
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
    rb, isbalanced = combo(rdr, cA, cB, bulk, args)

    if args['gccorr']:
        log('Applying additional GC correction')
        rdr = addgccorrect(rb, isbalanced, args)

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
    isbalanced = defaultdict(lambda : dict())
    cells = set(e for c in rdr for b in rdr[c] for e in rdr[c][b])
    jobs = ((c, b, np.random.randint(1000)) for c in rdr for b in rdr[c])
    njobs = sum(len(rdr[c].keys()) for c in rdr)
    snps = {c : sorted(cA[c].keys()) if c in cA else [] for c in rdr}
    bar = ProgressBar(total=njobs, length=40, verbose=False)
    initargs = (snps, cA, cB, bulk, cells, args)
    pool = Pool(processes=min(args['j'], njobs), initializer=init, initargs=initargs)
    counts = (lambda c, b, e : rdr[c][b][e].items())
    for c, b, A, B, isbal in pool.imap_unordered(combine, jobs):
        rb[c][b] = {e : dict(counts(c, b, e) + [('BAF', (A[e], B[e]) if e in A else (0, 0))]) for e in cells}
        isbalanced[c][b] = isbal
        bar.progress(advance=True, msg="Combined bin {}:{}-{}".format(c, b[0], b[1]))
    pool.close()
    pool.join()
    return rb, isbalanced
    

def init(_snps, _cA, _cB, _bulk, _cells, args):
    global snps, cA, cB, bulk, cells, blocksize, restarts, boot, alpha, alpha_correct, maxerror, minerror, missingsnps, minimumsnps, phasecorr
    snps = _snps
    cA = _cA
    cB = _cB
    bulk = _bulk
    cells = _cells
    blocksize = args['blocksize']
    restarts = args['restarts']
    boot = args['bootstrap']
    alpha = args['significance']
    alpha_correct = args['significance'] / (20 * float(sum(len(snps[c]) for c in snps)))
    maxerror = args['maxerror']
    minerror = args['minerror']
    minimumsnps = args['minimumsnps']
    missingsnps = args['missingsnps']
    phasecorr = args['phasecorr']


def combine(job):
    c, b, seed = job
    np.random.seed(seed)
    L = bisect.bisect_left(snps[c], b[0])
    R = bisect.bisect_right(snps[c], b[1])
    snpdensity = 1000.0 * max(0.0, (float(R - L) / float(b[1] - b[0])) if b[1] > b[0] else 0.0)

    if L >= R or snpdensity < minimumsnps:
        if missingsnps is None:
            return c, b, dict(), dict(), False
        else:
            return c, b, {e : missingsnps[0] for e in cells}, {e : missingsnps[1] for e in cells}, False

    snpsLR = snps[c][L:R]
    assert all(snpsLR[x - 1] < snpsLR[x] if x > 0 else True for x, o in enumerate(snpsLR))
    
    que = deque(snpsLR)
    assert sorted(snpsLR) == list(que) and b[0] <= que[0] and que[-1] <= b[1]
    omap = {}
    for bk in range(b[0], b[1]+1, blocksize):
        while que and bk <= que[0] < bk + blocksize:
            o = que.popleft()
            omap[o] = bk
    assert set(snpsLR) == set(omap.keys())

    # Fix potential phasing errors within each block
    def count_block(block):
        if phasecorr:
            if len(block) == 1:
                return (bulk[c][block[0]][0], bulk[c][block[0]][1])
            assert all(p < q for p, q in zip(block[:-1], block[1:])), "Positions in block {} are not sorted!".format(block)
            backward = {p : (sum(bulk[c][q][0] for q in block[:x+1]), 
                                sum(bulk[c][q][1] for q in block[:x+1])) for x, p in enumerate(block[:-1])}
            forward = {p : (sum(bulk[c][q][0] for q in block[x+1:]), 
                            sum(bulk[c][q][1] for q in block[x+1:])) for x, p in enumerate(block[:-1])}
            dotest = (lambda p, flip : sm.stats.proportions_ztest(count=(backward[p][0], forward[p][0 if not flip else 1]),
                                                                nobs=(sum(backward[p]), sum(forward[p])))
                                            if 0 < (backward[p][0] + forward[p][0 if not flip else 1]) < (sum(backward[p]) + sum(forward[p]))
                                            else (0.0, 1.0))
            tstat, test1 = zip(*(dotest(p, False) for p in block[:-1]))
            _, test2 = zip(*(dotest(p, True) for p in block[:-1]))
            pvals = sorted(set(test1))
            evaluate = (lambda test1, test2 : False if test1 >= alpha_correct or test1 >= test2 else True)
            flips = [evaluate(*t) for t in zip(test1, test2)]
            if sum(flips) > 0:
                sel, _ = max(((x, abs(s)) for x, s in enumerate(tstat) if flips[x]), key=(lambda p : p[1]))
            else:
                sel = len(block) + 1
            dofix = {p : x > sel for x, p in enumerate(block)}
        else:
            dofix = {p : False for p in block}
        return (sum(bulk[c][p][0 if not dofix[p] else 1] for p in block),
                sum(bulk[c][p][1 if not dofix[p] else 0] for p in block))

    blocks = {bk : count_block(sorted(o for o in omap if omap[o] == bk)) for bk in set(omap.values())}
    blocks = {bk : blocks[bk] for bk in blocks if sum(blocks[bk]) > 0}

    allblocks = blocks.values()
    nis = [sum(block) for block in allblocks]
    xis = [block[0] if np.random.random() < 0.5 else block[1] for block in allblocks]
    beta = max((EM(ns=nis, xs=xis, start=np.random.randint(low=1, high=49)/100.0) for r in xrange(restarts)), key=(lambda x : x[1]))[0]
    if maxerror is None:
        thres = max(minerror, est_error(ns=nis, significance=alpha, restarts=restarts, bootstrap=boot))
    else:
        thres = maxerror

    if thres <= abs(beta - 0.5):
        minel = (lambda p : 0 if p[0] < p[1] else 1)
        sumpairs = (lambda I : reduce((lambda x, y : (x[0] + y[0], x[1] + y[1])), I))
        if minel(sumpairs(allblocks)) == 0:
            swap = {o : False if blocks[omap[o]][0] < blocks[omap[o]][1] else True for o in snpsLR}
        else:
            swap = {o : False if blocks[omap[o]][1] < blocks[omap[o]][0] else True for o in snpsLR}
        isbal = False
    else:
        bkswap = {bk : False if np.random.random() < 0.5 else True for bk in blocks}
        swap = {o : True if bkswap[omap[o]] else False for o in snpsLR}
        isbal = True
                
    A = reduce(inupdate, (Counter(cA[c][o] if not swap[o] else cB[c][o]) for o in snpsLR))
    B = reduce(inupdate, (Counter(cB[c][o] if not swap[o] else cA[c][o]) for o in snpsLR))
    return c, b, A, B, isbal


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


def addgccorrect(data, isbalanced, args):
    gcs = gccount(args['gccorr'], data, args['j'])
    gccorrect(data, gcs, isbalanced, args['j'])

    
def gccount(ref, data, j=56):
    jobs = set(data.keys())
    bar = ProgressBar(total=len(jobs), length=min(len(jobs), 40), verbose=False)
    initargs = (ref, data)
    pool = Pool(processes=min(j, len(jobs)), initializer=init_gccounting, initargs=initargs)
    progress = (lambda c : bar.progress(advance=True, msg="Counted chromosome {}".format(c)))
    gcresult = {c : {b : gcs[b] for b in gcs} for c, gcs in pool.imap_unordered(gccounting, jobs) if progress(c)}
    pool.close()
    pool.join()
    return gcresult
    

def init_gccounting(_ref, _bins):
    global ref, bins
    ref = _ref
    bins = _bins


def gccounting(c):
    chro = c.replace('chr', '')
    sel = sorted(bins[c].keys(), key=(lambda b : (b[0], b[1])))
    tgts = deque(sel + [None])
    getchr = (lambda p : p[0].replace('>', '').replace('chr', ''))
    start = False
    tgt = tgts.popleft()
    offset = 0
    res = {b : Counter() for b in sel}
    with open(ref, 'r') as i:
        for l in i:
            if not start and l[0] == '>' and getchr(l.strip().split()) == chro:
                start = True
            elif start and l[0] == '>':
                break
            elif start:
                read = l.strip()
                while tgt is not None:
                    if tgt[1] <= offset + len(read):
                        if offset < tgt[1]:
                            res[tgt].update(read[max(0, tgt[0] - offset):tgt[1] - offset])
                        tgt = tgts.popleft()
                    else:
                        if tgt[0] < offset + len(read):
                            res[tgt].update(read[max(0, tgt[0] - offset):])
                        break
                offset += len(read)
    return c, res


def gccorrect(data, gcs, isbalanced, j):
    jobs = set(e for c in data for b in data[c] for e in data[c][b])
    bar = ProgressBar(total=len(jobs), length=min(len(jobs), 40), verbose=False)
    initargs = (data, gcs, isbalanced)
    pool = Pool(processes=min(j, len(jobs)), initializer=init_gccorrecting, initargs=initargs)
    progress = (lambda e : bar.progress(advance=True, msg="GC Corrected cell {}".format(e)))
    for e, gcs in pool.imap_unordered(gccorrecting, jobs):
        progress(e)
        for c in gcs:
            for b in gcs[c]:
                data[c][b][e]['RDR'] = gcs[c][b]
    pool.close()
    pool.join()
    return


def init_gccorrecting(_data, _gcs, _isbalanced):
    global data, gcs, isbalanced
    data = _data
    gcs = _gcs
    isbalanced = _isbalanced


def gccorrecting(e):
    getgc = (lambda C : ((C['C'] + C['G']) / float(C['C'] + C['G'] + C['A'] + C['T'])) if (C['C'] + C['G'] + C['A'] + C['T']) > 0 else 0.5)
    form = (lambda c, b : {'Tumor' : data[c][b][e]['readcount'], 'Normal' : data[c][b][e]['normalcount'], '%GC' : getgc(gcs[c][b]), 'isbalanced' : isbalanced[c][b] if isbalanced is not None else True})
    curr = {(c, b[0], b[1]) : form(c, b) for c in data for b in data[c]}
    with warnings.catch_warnings() as w:
        warnings.simplefilter("ignore")
        tcorr = gccorr(curr, rkey='Tumor')
        ncorr = gccorr(curr, rkey='Normal')
    getrdr = (lambda c, b : (tcorr[(c, b[0], b[1])] / ncorr[(c, b[0], b[1])]) if ncorr[(c, b[0], b[1])] > 0 else 0.0)
    return e, {c : {b : getrdr(c, b) for b in data[c]} for c in data}


def gccorr(curr, rkey='Tumor', plot=False, frac=0.3, tol=0.1, genomethres=0.4):
    reg = (lambda D : D / (np.sum(D) / float(len(D))))
    regD = (lambda D : {b : D[b] / (sum(D.values()) / float(len(D))) for b in D})
    ixs = sorted(curr.keys(), key=(lambda b : curr[b]['%GC']))

    sel1 = [b for b in ixs if curr[b]['isbalanced']]
    if (len(sel1) / float(len(ixs))) <= genomethres:
        sel1 = [b for b in ixs]
    thres = int(round(len(sel1) * 0.0))
    firlas = (lambda L : (L[0], L[-1]))
    min1, max1 = firlas(sorted([curr[b][rkey] for b in sel1])[thres:len(sel1)-thres])
    xs1, ys1 = map(np.array, zip(*((curr[b]['%GC'], curr[b][rkey]) for b in sel1 if min1 <= curr[b][rkey] <= max1)))
    ys1 = reg(ys1)
    id1 = ys1 >= tol
    if not any(id1):
        return regD({b : curr[b][rkey] for b in ixs})
    xs1, ys1 = xs1[id1], ys1[id1]
    lowxs1, lowys1 = sm.nonparametric.lowess(ys1, xs1, frac=frac).T
    fs1 = np.interp(x=xs1, xp=lowxs1, fp=lowys1)
    if plot:
        plt.scatter(xs1, ys1, c='cyan', s=20)
        plt.plot(xs1, fs1, c='blue')

    ref, mid = max(zip(xs1, fs1), key=(lambda p : p[1]))
    if np.isnan(mid):
        return regD({b : curr[b][rkey] for b in ixs})
    mar = 2.0 * np.max(fs1)
    mir = np.min(fs1) / 2.0
    assert mir <= mid <= mar, (mir, mid, mar)
    dist = (lambda s : np.count_nonzero(np.abs(ys1 - (fs1 + s)) <= tol))
    fs1 = fs1 + max(np.linspace(mir, mar + 1, 100) - mid, key=dist)
    id2 = np.abs(fs1 - ys1) <= tol
    if not any(id2):
        return regD({b : curr[b][rkey] for b in ixs})
    xs2, ys2 = xs1[id2], ys1[id2]
    lowxs2, lowys2 = sm.nonparametric.lowess(ys2, xs2, frac=frac).T
    fs2 = np.interp(x=xs2, xp=lowxs2, fp=lowys2)
    if plot:
        plt.scatter(xs2, ys2, c='yellow', s=20)
        plt.plot(xs2, fs2, c='orange')

    fs = np.interp(x=xs1, xp=xs2, fp=fs2)
    xs = map((lambda b : curr[b]['%GC']), ixs)
    assert xs == sorted(xs), xs
    fs = reg(np.interp(x=xs, xp=xs1, fp=fs))
    fs[(np.isnan(fs)) | (fs < tol)] = tol
    res = regD({b : curr[b][rkey] * (max(fs) / f) for b, f in zip(ixs, fs)})
    if plot:
        plt.scatter(xs, reg(map((lambda b : curr[b][rkey]), ixs)), c='k', s=10)
        plt.plot(xs, fs, c='green')
        plt.scatter(xs, map((lambda b : res[b]), ixs), c='red', s=10)
    
    return res
    

    
if __name__ == '__main__':
    main()

