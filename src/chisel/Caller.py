#!/usr/bin/env python2.7

import sys, os
import argparse
import math
import itertools
import multiprocessing as mp

from multiprocessing import Lock, Value, Pool
from collections import defaultdict
from collections import Counter

import numpy as np
import scipy.stats

from Utils import *
from Clusterizer import *
from Combiner import EM


def parse_args(args):
    description = "Inferring allele- and haplotype-specific copy number in single cells from estimated RDRs and BAFs."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", type=str, help="Input file with combined RDR and BAF per bin and per cell")
    parser.add_argument("-A","--sensitivity", required=False, type=float, default=1.0, help="Sensitivity of model selection for ploidy (default: 1, increase this parameter to lower sensitivity to noisy data, adjust this value (e.g. 2, 4, ..., 10, ...) to better deal with high-variance data (e.g. low coverage, small number of cells, low number of phased SNPs, etc...)")
    parser.add_argument("-K","--upperk", required=False, type=int, default=100, help="Maximum number of bin clusters (default: 100, use 0 to consider maximum number of clusters)")
    parser.add_argument("-P","--maxploidy", required=False, type=int, default=4, help="Maximum total copy number to consider for base cluster (default: 4, corresponding to 2 WGD)")
    parser.add_argument("-k","--lowerk", required=False, type=int, default=0, help="Minimum number of bin clusters (default: 0)")
    parser.add_argument("-t","--threshold", required=False, type=float, default=0.05, help="Maximum tolerated difference of unexplained variance for cluster selection (default: 0.05)")
    parser.add_argument("-q","--restarts", required=False, type=int, default=200, help="Number of restarts for clustering (default: 200)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-d","--shift", required=False, type=float, default=0.05, help="Maximum estimate shift for base cluster, the value is increased when base is not found (default: 0.05)")
    parser.add_argument("-s","--significativity", required=False, type=float, default=0.02, help="Minimum proportion of the clusters use to select ploidy (default: 0.02)")
    parser.add_argument("-l","--lorder", required=False, type=int, default=1, help="Order of the l-norm distance to be used, either 1 or 2 (default: 1)")
    parser.add_argument("--fastscaling", required=False, default=False, action='store_true', help="Consider average BAF of clusters instead of EM (default: False, using it is generally safe and significantly increases speed)")
    #parser.add_argument("--scoring", required=False, default=False, action='store_true', help=" (default: equal to number of available processors)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: none)")
    args = parser.parse_args(args)

    if not os.path.isfile(args.INPUT):
        raise ValueError("Input file does not exist: {}".format(args.INPUT))

    if not 0.0 <= args.significativity <= 1.0:
        raise ValueError("Significativity must be in [0, 1]!")
    if args.maxploidy < 2:
        raise ValueError("Maximum total copy number for base must be at least 2!")
    if not 0.0 <= args.shift <= 0.5:
        raise ValueError("Base shift must be in [0, 0.5]!")
    if args.upperk < 0:
        raise ValueError("Maximum number of clusters must be positive or zero!")
    if args.lowerk < 0:
        raise ValueError("Minimum number of clusters must be positive or zero!")
    if not 0.0 <= args.threshold <= 1.0:
        raise ValueError("Threshold must be in [0, 1]!")
    if args.lorder not in {1, 2}:
        raise ValueError("L-norm order must be either 1 or 2!")
    if args.restarts <= 0:
        raise ValueError("restarts must be positive!")
    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")
    if args.seed and args.seed < 0:
        raise ValueError("Random seed must be positive or zero!")
    else:
        np.random.seed(args.seed)

    return {
        'input' : args.INPUT,
        'sensitivity' : args.sensitivity,
        'significativity' : args.significativity,
        'maxploidy' : args.maxploidy,
        'shift' : args.shift,
        'UB' : args.upperk,
        'LB' : args.lowerk,
        'e' : args.threshold,
        'lord' : args.lorder,
        'restarts' : args.restarts,
        'j' : args.jobs,
        'fastscaling' : args.fastscaling,
        #'scoring' : args.scoring,
        'scoring' : False,
        'seed' : args.seed
    }


def main(args=None, stdout_file=None):
    log('Parsing and checking arguments')
    args = parse_args(args)
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]), level='INFO')

    log('Reading combined RDRs and BAFs of barcoded cells')
    bins = read_combo(args['input'])

    log('Formatting RDRs and BAFs')
    data, mapb, cells = format_combo(bins)

    log('Clustering globally')
    par = {}
    par['data'] = data
    par['restarts'] = args['restarts']
    par['threshold'] = args['e']
    par['seed'] = args['seed']
    par['lord'] = args['lord']
    par['j'] = args['j']
    par['LB'] = args['LB']
    par['UB'] = args['UB']
    members = kclustering(**par)

    log('Estimating RDR and BAF of every cluster')
    clusters = est_clusters(bins, cells, members, mapb)
    
    log('Selecting ploidies')
    calls, ploidy = call(bins, members, mapb, clusters, cells, args)
    msg = '\n'.join(['Cells with base ploidy {}: {}'.format(p, v) for p, v in Counter(ploidy.values()).items()])
    log('Number of cells for every ploidy\' level:\n{}'.format(msg), level='INFO')

    log('Inferring copy numbers')
    copynumbers = infer(bins, cells, mapb, members, calls)

    log('Phasing copy-number states along the genome')
    copynumbers = phase(copynumbers, cells, args['j'])

    log('Writing results')
    order = (lambda b : (orderchrs(b[0]), int(b[1]), int(b[2])))
    form_bins = (lambda b, e : [bins[b][e]['NORM'], bins[b][e]['wRDR'], bins[b][e]['RDR'], bins[b][e]['AB'][0], bins[b][e]['AB'][1], bins[b][e]['BAF']])
    form_clus = (lambda b, e : [int(members[mapb[b]]), '{}|{}'.format(*copynumbers[b][e])])
    form = (lambda b, e : '\t'.join(map(str, [b[0], b[1], b[2], e] + form_bins(b, e) + form_clus(b, e))))
    gen = (form(b, e) for b in sorted(bins, key=order) for e in sorted(cells))
    header = '\t'.join(['#CHR', 'START', 'END', 'CELL', 'NORM_COUNT', 'COUNT', 'RDR', 'A_COUNT', 'B_COUNT', 'BAF', 'CLUSTER', 'CN_STATE'])

    if stdout_file is not None:
        with open(stdout_file, 'w') as f:
            f.write(header + '\n')
            f.write('\n'.join(gen))
    else:
        print header
        for g in gen:
            print g


def read_combo(f):
    form = (lambda p : ((p[0], int(p[1]), int(p[2])), p[3], int(p[4]), int(p[5]), float(p[6]), int(p[7]), int(p[8]), float(p[9])))
    data = defaultdict(lambda : dict())
    with open(f, 'r') as i:
        for l in (l for l in i if l[0] != '#' and len(l) > 1):
            b, e, norm, wrdr, rdr, A, B, baf = form(l.strip().split())
            assert e not in data[b], 'A cell with duplicated name has been found for the same bin {}'.format(b)
            data[b][e] = {'RDR':rdr, 'NORM':norm, 'wRDR':wrdr, 'BAF':baf, 'AB':(A, B), 'wBAF':A+B, 'mirroredBAF':min(baf, 1.0-baf)}
    return {b : data[b] for b in data}


def format_combo(bins):
    cells = set(e for b in bins for e in bins[b])
    if False in set(cells == set(bins[b].keys()) for b in bins):
        raise ValueError('ERROR: RDR and BAF for some cells are not defined across the same set of bins!')
    cells = sorted(cells)

    order = (lambda b : (orderchrs(b[0]), int(b[1]), int(b[2])))
    pos = sorted(bins.keys(), key=order)
    mapb = {b : x for x, b in enumerate(pos)}
    data = np.array([[bins[b][e][d] for e in cells for d in ('RDR', 'mirroredBAF')] for b in pos])

    return data, mapb, cells


def est_clusters(bins, cells, members, mapb):
    partition = {c : filter((lambda b : members[mapb[b]] == c), bins.keys()) for c in set(members)}
    assert len(members) == sum(len(partition[c]) for c in partition), 'The bins have not been partitioned into clusters, and something is missing'
    partition = {c : partition[c] for c in partition if len(partition[c]) > 0}
    size = {c : sum(b[2] - b[1] for b in partition[c]) for c in partition}
    cal_rdr = (lambda c, e : sum(bins[b][e]['RDR'] * bins[b][e]['wRDR'] for b in partition[c]) / float(sum(bins[b][e]['wRDR'] for b in partition[c])))
    est_rdr = (lambda c, e : cal_rdr(c, e) if sum(bins[b][e]['wRDR'] for b in partition[c]) > 0 else 0.0)
    cal_baf = (lambda c, e : sum(bins[b][e]['mirroredBAF'] * bins[b][e]['wBAF'] for b in partition[c]) / float(sum(bins[b][e]['wBAF'] for b in partition[c])))
    est_baf = (lambda c, e : cal_baf(c, e) if sum(bins[b][e]['wBAF'] for b in partition[c]) > 0 else 0.5)
    est = (lambda c, e : {'RDR' : est_rdr(c, e), 'BAF' : est_baf(c, e), 'size' : size[c]})
    return {c : {e : est(c, e) for e in cells} for c in partition}


def call(bins, members, mapb, clusters, cells, args):
    jobs = (e for e in cells)
    bar = ProgressBar(total=len(cells), length=40, verbose=False)
    totsize = set(sum(clusters[c][e]['size'] for c in clusters) for e in cells)
    assert len(totsize) == 1, 'The same cluster has different sizes in different cells'
    counts = {e : {c : {tuple(b) : tuple(bins[b][e]['AB']) for b in bins if members[mapb[b]] == c} for c in clusters} for e in cells}
    rdrs = {e : {c : {tuple(b) : (bins[b][e]['RDR'], bins[b][e]['wRDR']) for b in bins if members[mapb[b]] == c} for c in clusters} for e in cells}
    # initargs = (cells, bins, members, mapb, clusters, list(totsize)[0], args['maxploidy'], args['significativity'], args['shift'], args['lord'])
    initargs = (clusters, list(totsize)[0], args, counts, rdrs)
    pool = Pool(processes=min(args['j'], len(cells)), initializer=init, initargs=initargs)
    calls = {e : None for e in cells}
    ploidy = {e : 0 for e in cells}
    progress = (lambda e : bar.progress(advance=True, msg="Called cell {}".format(e)))
    update = (lambda (e, cns, p) : not calls.update({e : cns}) and not ploidy.update({e : p}) and progress(e))
    if not args['scoring']:
        res = map(update, pool.imap_unordered(calling, jobs))
    else:
        res = map(update, pool.imap_unordered(scoring, jobs))
    assert False not in res, 'Some process in the pool did not terminate succesfully'
    pool.close()
    pool.join()
    return {c : {e : calls[e][c] for e in cells} for c in clusters}, ploidy


def init(_clusters, _totsize, args, _counts, _rdrs):
    global clusters, totsize, maxploidy, significativity, shift, lord, fastscaling, counts, rdrs, sensitivity
    clusters = _clusters
    totsize = _totsize
    maxploidy = args['maxploidy']
    sensitivity = args['sensitivity']
    significativity = args['significativity']
    shift = args['shift']
    lord = args['lord']
    fastscaling = args['fastscaling']
    counts = _counts
    rdrs = _rdrs


def scoring(e):
    pre_states = None
    pre_scores = None
    pre_ploidy = None
    ccnts = {c : {tuple(b) : tuple(counts[e][c][b]) for b in counts[e][c]} for c in counts[e]}

    for ploidy in range(2, maxploidy+1, 2):
        scale = scaling(e, ploidy, clusters, ccnts, shift)
        ctot = {c : int(round(clusters[c][e]['RDR'] * scale)) for c in clusters}
        cmax = max(ctot.values())
        split = (lambda c : (ctot[c] * (1.0 - clusters[c][e]['BAF']), ctot[c] * clusters[c][e]['BAF']))
        clus = {c : split(c) for c in clusters}
        assert False not in set(clus[c][0] >= clus[c][1] for c in clusters), 'Found wrong order in allele-specific copy numbers for scoring'

        mkstate = (lambda t, s : (max(s, t - s), min(s, t - s)))
        allstates = {t : set(mkstate(t, s) for s in xrange(t+1)) for t in xrange(cmax+1)}
        assert False not in set(sum(s) <= t and s[0] >= s[1] for t in allstates for s in allstates[t]), 'Found a state with maximum copy number higher than predicted maximum'

        dist = (lambda x, y : np.linalg.norm(np.array(x) - np.array(y), ord=lord))
        states = {c : min(allstates[ctot[c]], key=(lambda s : dist(clus[c], s))) for c in clusters}
        scores = {s : sum(clusters[c][e]['size'] for c in filter((lambda c : states[c] == s), states.keys())) for s in set(states.values())}

        assert sum(clusters[c][e]['size'] for c in states) == sum(scores.values()), 'The scoring value did not account for all clusters according to their size'
        scores = {s : scores[s] for s in scores if (float(scores[s]) / float(totsize)) >= significativity}
        if pre_scores and len(pre_scores.keys()) >= len(scores.keys()):
            return e, pre_states, pre_ploidy
        else:
            pre_states = states
            pre_scores = scores
            pre_ploidy = ploidy

    assert states == pre_states, 'The states with the highest ploidy have not been considered'
    return e, pre_states, pre_ploidy


def calling(e):
    ccnts = {c : {tuple(b) : tuple(counts[e][c][b]) for b in counts[e][c]} for c in counts[e]}
    crdrs = {c : {tuple(b) : tuple(rdrs[e][c][b]) for b in rdrs[e][c]} for c in rdrs[e]}
    assert False not in set(set(ccnts[c].keys()) == set(crdrs[c].keys()) for c in ccnts), 'RDRs or BAFs for some clusters is missing'
    baf = (lambda A, B : (float(min(A, B)) / float(A + B)) if (A + B) > 0 else 0.5)
    cbafs = {c : {b : baf(*ccnts[c][b]) for b in ccnts[c]} for c in ccnts}

    rdrbase, selected = identify_base(e, clusters, ccnts, shift)
    scale = {pl : float(pl) / rdrbase for pl in range(2, maxploidy+1, 2)}

    norm = (lambda S : float((len(S) - 1) if len(S) > 1 else 1))
    mean = (lambda L : float(sum(e for e in L)) / float(len(L)))
    suvar = (lambda L, u : float(sum(math.pow(e - u, 2) for e in L)) / norm(L))
    var = (lambda L : suvar(L, mean(L)) if len(L) > 1 else 0.001)
    var = (lambda L, u : suvar(L, u) if len(L) > 0 else 0.001)

    avg = (lambda W : float(sum(w for w in W)) / float(len(W)))
    cmax = {pl : max(int(round(clusters[c][e]['RDR'] * scale[pl])) for c in clusters) for pl in scale}
    mkstate = (lambda t, s : (max(s, t - s), min(s, t - s)))
    allstates = {pl : {t : set(mkstate(t, s) for s in xrange(1, t+1)) for t in xrange(1, cmax[pl]+1)} for pl in scale}
    apxtot = {pl : {c : max(1, int(round(clusters[c][e]['RDR'] * scale[pl]))) for c in clusters} for pl in scale}
    apxsta = {pl : {c : min(allstates[pl][apxtot[pl][c]], key=(lambda s : abs(baf(*s) - clusters[c][e]['BAF']))) for c in clusters} for pl in scale}
    flatst = {pl : set(s for t in allstates[pl] for s in allstates[pl][t]) for pl in scale}
    bafvar = {pl : {s : var([cbafs[c][b] for c in clusters if apxsta[pl][c] == s for b in cbafs[c]], baf(*s)) for s in flatst[pl]} for pl in scale}
    rdrvar = {pl : {s : var([crdrs[c][b][0] for c in clusters if apxsta[pl][c] == s for b in crdrs[c]], sum(s)/scale[pl]) for s in flatst[pl]} for pl in scale}
    mbafvar = max(avg([bafvar[pl][s] for s in bafvar[pl]]) for pl in scale)
    mrdrvar = max(avg([rdrvar[pl][s] for s in rdrvar[pl]]) for pl in scale)

    gauss = scipy.stats.norm.logpdf

    baflh = (lambda x, var, A, B : gauss(x=x, loc=baf(A, B), scale=math.sqrt(var)))
    rdrlh = (lambda x, var, A, B, pl : gauss(x=x, loc=(float(A+B)/scale[pl]), scale=math.sqrt(var)))
    blh = (lambda b, vb, r, vr, A, B, pl : baflh(b, vb, A, B) + rdrlh(r, vr, A, B, pl))
    prep = (lambda A, B, c, pl, vb, vr : sum(blh(cbafs[c][b], vb, crdrs[c][b][0], vr, A, B, pl) for b in ccnts[c]))
    vbaf = (lambda u, c, pl : max(mbafvar, sum(math.pow(cbafs[c][b] - u, 2) for b in cbafs[c]) / norm(cbafs[c]) ))
    vrdr = (lambda u, c, pl : max(mrdrvar, sum(math.pow(crdrs[c][b][0] - u, 2) for b in crdrs[c]) / norm(crdrs[c])))
    clh = (lambda A, B, c, pl : prep(A, B, c, pl, vbaf(baf(A, B), c, pl), vrdr(float(A+B)/scale[pl], c, pl)))

    cmax = {pl : max(int(round(clusters[c][e]['RDR'] * scale[pl])) for c in clusters) for pl in scale}
    mkstate = (lambda t, s : (max(s, t - s), min(s, t - s)))
    allstates = {pl : {t : set(mkstate(t, s) for s in xrange(1, t+1)) for t in xrange(1, cmax[pl]+1)} for pl in scale}
    capp = {pl : {c : int(round(clusters[c][e]['RDR'] * scale[pl])) for c in clusters} for pl in scale}
    bnd = (lambda v, pl : max(1, min(cmax[pl], v)))
    avail = (lambda pl, c : set(st for t in set(bnd(capp[pl][c]+i*d, pl) for i in xrange(2) for d in [1, -1]) for st in allstates[pl][t]))
    vmax = (lambda L, c, pl : max(((l, clh(l[0], l[1], c, pl)) for l in L), key=(lambda x : x[1])))
    states = {pl : {c : vmax(avail(pl, c), c, pl) if clusters[c][e]['RDR'] > 0 else ((0, 0), 1.0) for c in clusters} for pl in scale}
    
    norm = float(sum(clusters[c][e]['size'] for c in clusters if clusters[c][e]['RDR'] > 0))
    issupp = (lambda pl, t, THRES : (sum(clusters[c][e]['size'] for c in clusters if t == sum(states[pl][c][0])) / norm) >= THRES)
    safemax = (lambda L : max(L) if len(L) > 0 else max(t for t in allstates[pl] if issupp(pl, t, 0.00)))
    maxsupp = {pl : safemax([t for t in allstates[pl] if issupp(pl, t, 0.02)]) for pl in states}
    n = sum(len(ccnts[c]) for c in ccnts)
    lh = {pl : sum(states[pl][c][1] for c in states[pl]) for pl in states}
    ps = {pl : 1 + sum(len(allstates[pl][t]) for t in allstates[pl] if t <= maxsupp[pl]) for pl in states}
    bic = (lambda pl : sensitivity * math.log(n) * float(ps[pl]) * 2 - 2.0 * lh[pl])
    best = min(lh.keys(), key=bic)

    return e, {c : states[best][c][0] for c in clusters}, best


def identify_base(e, clusters, ccnts, shift):
    nis = (lambda c : [sum(ccnts[c][b]) for b in ccnts[c]])
    xis = (lambda c : [ccnts[c][b][0] if np.random.random() < 0.5 else ccnts[c][b][1] for b in ccnts[c]])
    gstart = (lambda : np.random.randint(low=1, high=49) / 100.0)
    beta = (lambda c, nis, xis : max((EM(ns=nis, xs=xis, start=gstart()) for r in xrange(50)), key=(lambda x : x[1]))[0])
    wavg = (lambda L : sum(l[0] * l[1] for l in L) / float(sum(l[1] for l in L)))
    s = shift
    while s <= 0.5:
        if not fastscaling:
            isneutral = (lambda c : abs(0.5 - clusters[c][e]['BAF']) <= (0.2 + s) and clusters[c][e]['RDR'] > 0.0 and abs(0.5 - beta(c, nis(c), xis(c))) <= s)
        else:
            isneutral = (lambda c : abs(0.5 - clusters[c][e]['BAF']) <= s and clusters[c][e]['RDR'] > 0.0)
        sel = filter(isneutral, clusters.keys())
        if len(sel) > 0:
            rec = (lambda c : (clusters[c][e]['RDR'], clusters[c][e]['size']))
            mk_bubble = (lambda v : set(rec(c) for c in sel if abs(v - clusters[c][e]['RDR']) <= 0.05))
            collapse = (lambda bu : (wavg(bu), sum(u[1] for u in bu)))
            return max((collapse(mk_bubble(clusters[c][e]['RDR'])) for c in sel), key=(lambda x : x[1]))[0], sel
        else:
            s += shift/2.0
    assert False, 'There is a bin with a BAF shift > 0.5, likely BAF was not mirrored between 0 and 0.5'


def infer(bins, cells, mapb, members, calls):
    reg = (lambda b, e, p : (max(p), min(p)) if bins[b][e]['BAF'] <= 0.5 else (min(p), max(p)))
    getcn = (lambda b, e : reg(b, e, calls[members[mapb[b]]][e]))
    return {b : {e : getcn(b, e) for e in cells} for b in bins}


def phase(cns, cells, j):
    jobs = set(b[0] for b in cns)
    bar = ProgressBar(total=len(jobs), length=40, verbose=False)
    order = (lambda b : (orderchrs(b[0]), int(b[1]), int(b[2])))
    initargs = ({c : {b : {e : cns[b][e] for e in cells} for b in cns if b[0] == c} for c in jobs}, sorted(cells))
    pool = Pool(processes=min(j, len(jobs)), initializer=init_phasing, initargs=initargs)
    progress = (lambda c : bar.progress(advance=True, msg="Phased chromosome {}".format(c)))
    return {b : phased[b] for c, phased in pool.imap_unordered(phasing, jobs) if progress(c) for b in phased}


def init_phasing(_table, _cells):
    global table, cells
    table = _table
    cells = _cells


def phasing(c):
    D = []
    B = []
    order = (lambda b : (orderchrs(b[0]), int(b[1]), int(b[2])))
    swapped = (lambda p : (p[1], p[0]))
    norm = {i : {e : table[c][b][e] for e in cells} for i, b in enumerate(sorted(table[c].keys(), key=order))}
    swap = {i : {e : swapped(table[c][b][e]) for e in cells} for i, b in enumerate(sorted(table[c].keys(), key=order))}
    hd = (lambda x, y : 1 if x != y else 0)
    dist = (lambda x, y : sum(hd(x[e][0], y[e][0]) + hd(x[e][1], y[e][1]) for e in cells))

    for i in xrange(len(table[c].keys())):
        if i == 0:
            Dnon, Dswa = 0, 0
            Bnon, Bswa = -1, -1
        else:
            Dnonnon = D[i - 1]['0'] + dist(norm[i - 1], norm[i])
            Dswanon = D[i - 1]['1'] + dist(swap[i - 1], norm[i])
            Dnon = min(Dnonnon, Dswanon)
            Bnon = '0' if Dnonnon <= Dswanon else '1'

            Dnonswa = D[i - 1]['0'] + dist(norm[i - 1], swap[i])
            Dswaswa = D[i - 1]['1'] + dist(swap[i - 1], swap[i])
            Dswa = min(Dnonswa, Dswaswa)
            Bswa = '0' if Dnonswa <= Dswaswa else '1'

        D.append({'0' : Dnon, '1' : Dswa})
        B.append({'0' : Bnon, '1' : Bswa})

    res = []
    for i in reversed(xrange(len(table[c]))):
        if i == len(table[c]) - 1:
            if D[i]['0'] < D[i]['1']:
                res.insert(0, norm[i])
                prev = B[i]['0']
            else:
                res.insert(0, swap[i])
                prev = B[i]['1']
        else:
            res.insert(0, norm[i] if prev == '0' else swap[i])
            prev = B[i][prev]

    assert prev == -1, 'The backtracking of the dynamic programming did not terminate with the first bin'
    return c, {b : res[i] for i, b in enumerate(sorted(table[c].keys(), key=order))}


if __name__ == '__main__':
    main()
