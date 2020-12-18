#!/usr/bin/env python2.7

import sys, os
import argparse
import shutil
import random
import warnings

from itertools import cycle
from collections import defaultdict
from collections import Counter

import numpy as np
import scipy
import scipy.cluster
import scipy.cluster.hierarchy as hier

from Utils import *


def parse_args():
    description = "Infer clones as subpopulations of cells with the same complement of CNAs and outputs a file with the mapping of every cell to the corresponding clone."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", type=str, help="Input file with RDR, BAF, and inferred copy numbers.")
    parser.add_argument("-f", "--maxdiff", required=False, type=float, default=0.06, help="Maximum fraction of the genome with different copy-number states allowed in a clone (default: 0.06, when -1 is chosen the maximum cluster method of SciPy is used)")
    parser.add_argument("-r", "--refinement", required=False, type=float, default=0.0, help="Maximum difference to assign noisy cells to a clone (default: 0.0)")
    parser.add_argument("-s", "--minsize", required=False, type=int, default=14, help="Minimum size of subpopultation to define a clone (default: 14)")
    parser.add_argument("-l", "--linkage", required=False, type=str, default='single', help="Linkage method to use for the hierarchical clustering (default: single, it must be a valid linkage method available in SciPy when using a non-euclidean distance, i.e. 'single', 'complete', 'average', 'weighted')")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: none)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError('ERROR: input file does not exist!')
    if (not 0.0 <= args.maxdiff <= 1.0) and args.maxdiff != -1:
        raise ValueError('ERROR: the maximum different fraction of the genome must be either within [0, 1] or equal to -1!')
    if args.refinement is None:
        args.refinement = args.maxdiff
    if not 0.0 <= args.refinement <= 1.0:
        raise ValueError('ERROR: the refinement must be either within [0, 1]!')
    if args.minsize <= 0:
        raise ValueError('ERROR: the minimum size of subpopulations must be positive!')
    if not args.linkage in {'single', 'complete', 'average', 'weighted'}:
        raise ValueError('ERROR: the linkage method is invalid or not available for non-euclidean distances!')
    if args.seed and args.seed < 0:
        raise ValueError("Random seed must be positive or zero!")
    else:
        random.seed(args.seed)
        np.random.seed(args.seed)

    return {
        'input' : args.INPUT,
        'maxdiff' : args.maxdiff,
        'refinement' : args.refinement,
        'minsize' : args.minsize,
        'linkage' : args.linkage,
        'seed' : args.seed
    }


def main():
    log('Parsing and checking arguments')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]), level='INFO')

    log('Reading input')
    cns, pos, cells = reading(args['input'])

    log('Clustering cells in clones')
    clus = clustering(cns, pos, cells, args['maxdiff'], args['linkage'])

    log('Selecting clones')
    clones = selecting(clus, args['minsize'])
    log('Number of identified clones: {}'.format(len(set(clones.values()))), level='INFO')

    if len(clones) > 0 and args['refinement'] > 0.0:
        log('Refining clustering')
        clones, clus = refining(cns, clus, clones, args['refinement'])
    log('Number of discarded cells: {} over {} in total'.format(len(set(cells) - set(clones.keys())), len(set(cells))), level='INFO')

    log('Profiling clones')
    profiles = profiling(cns, clus)
    
    log('Writing clone map')
    print '\t'.join(['#CELL', 'CLUSTER', 'CLONE'])
    for c in cells:
        print '\t'.join(map(str, [c, clus[c], 'Clone{}'.format(clones[c]) if c in clones else 'None']))

    log('Writing clone-corrected copy numbers in provided input')
    ftmp = args['input'] + '_TMP'
    assert not os.path.isfile(ftmp), "Temporary file {} does already exist!".format(ftmp)
    form = (lambda p : ((p[0], int(p[1]), int(p[2])), p[3], p[0:12]))
    with open(args['input'], 'r') as i:
        with open(ftmp, 'w') as o:
            for l in i:
                if '#' != l[0]:
                    b, e, val = form(l.strip().split())
                    o.write('\t'.join(val + ['{}|{}'.format(*profiles[b][clus[e]])]) + '\n')
                else:
                    o.write('\t'.join(['#CHR', 'START', 'END', 'CELL', 'NORM_COUNT', 'COUNT', 'RDR', 'A_COUNT', 'B_COUNT', 'BAF', 'CLUSTER', 'HAP_CN', 'CORRECTED_HAP_CN']) + '\n')
    shutil.move(ftmp, args['input'])


def reading(f):
    cns = defaultdict(lambda : dict())
    form = (lambda p : ((p[0], int(p[1]), int(p[2])), p[3], tuple(map(int, p[11].split('|')))))
    with open(f, 'r') as i:
        for l in i:
            if l[0] != '#' and len(l) > 1:
                b, c, cn = form(l.strip().split())
                assert c not in cns[b] # and c not in stuff[b]
                cns[b][c] = cn
    cns = dict(cns)
    orderchrs = (lambda x : int(''.join([l for l in x if l.isdigit()])))
    order = (lambda b : (orderchrs(b[0]), int(b[1]), int(b[2])))
    pos = sorted(cns.keys(), key=order)
    cells = sorted(set(c for b in cns for c in cns[b]))
    return cns, pos, cells


def clustering(cns, pos, cells, maxdiff, linkage):
    states = {s : x for x, s in enumerate(set(cns[b][c] for b in pos for c in cells))}
    data = [[states[cns[b][c]] for b in pos] for c in cells]
    linkage = hier.linkage(data, method=linkage, metric='hamming', optimal_ordering=True)
    if maxdiff != -1:
        clus = hier.fcluster(linkage, t=maxdiff, criterion='distance')
    else:
        clus = hier.fcluster(linkage, t=len(cells), criterion='maxclust')
    return {e : clus[i] for i, e in enumerate(cells)}


def selecting(clus, minsize):
    size = {i : sum(clus[c] == i for c in clus) for i in set(clus.values())}
    return {c : clus[c] for c in clus if size[clus[c]] >= minsize}


def refining(cns, clus, chosen, maxdiff):
    clones = set(chosen.values())
    safeargmax = (lambda C : argmax(C) if len(C) > 0 else (1, 1))
    getcn = (lambda g, i : safeargmax(Counter([cns[g][c] for c in chosen if chosen[c] == i])))
    profile = {g : {i : getcn(g, i) for i in clones} for g in cns}
    diff = (lambda i, c, g : 1 if profile[g][i] != cns[g][c] else 0)
    weight = (lambda i, c : float(sum(diff(i, c, g) for g in profile)) / float(len(profile)))
    closest = (lambda c : min([(i, weight(i, c)) for i in clones], key=(lambda x : x[1])))
    ref = {c : closest(c) for c in clus if c not in chosen.keys()}
    newclones = {c : chosen[c] if c in chosen else ref[c][0] for c in clus if c in chosen or ref[c][1] <= maxdiff}
    newclus = {c : newclones[c] if c in newclones else clus[c] for c in clus}
    assert False not in set(len({clus[c], chosen[c], newclus[c], newclones[c]}) == 1 for c in chosen)
    return newclones, newclus


def profiling(cns, clus):
    clones = set(clus.values())
    mapclo = {i : filter(lambda e : clus[e] == i, clus.keys()) for i in clones}
    assert all(len(mapclo[i]) > 0 for i in mapclo), 'Found cluster assignment with no corresponding cell'
    getcn = (lambda g, i : argmax(Counter([cns[g][e] for e in mapclo[i]])))
    return {g : {i : getcn(g, i) for i in clones} for g in cns}    


if __name__ == '__main__':
    main()
