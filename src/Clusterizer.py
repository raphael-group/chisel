#!/usr/bin/env python2.7

import os, sys
import argparse
import math
import ctypes
import warnings

import multiprocessing as mp
from multiprocessing import Lock, Value, Pool

import numpy as np
import scipy.spatial

from Utils import *



def kclustering(data, restarts, threshold, seed=None, lord=1, j=1, LB=1, UB=None):
    pool_size = mp.cpu_count()
    os.system('taskset -cp 0-%d %s > /dev/null' % (pool_size, os.getpid()))

    error, center, pdist, cdist = getord(lord)
    if len(set(len(p) for p in data)) != 1:
        raise ValueError('All points must have the same length')
    if seed:
        np.random.seed(seed)
    
    points = np.array(data)
    TCENTER = center(points)
    TERROR = sum(error(p, TCENTER) for p in points)
    PAIRWISE = pdist(points)
    assert np.isfinite(PAIRWISE).all(), 'Pairwise distance contain NaN!\n{}'.format(PAIRWISE)

    objs = {}
    clus = {}
    
    if UB:
        R = min(len(data), UB)
    else:
        R = len(points)
        objs[R] = 0.0
        clus[R] = [i for i, p in enumerate(points)]

    if LB:
        L = max(0, LB)
    else:
        L = 0
        objs[L] = 1.0
        clus[L] = [0 for p in points]
    
    def compute(K):
        log('Computing for {}:'.format(K), level='INFO')
        if K not in objs:
            assert K not in clus, 'The number of clusters {} does not have an objective but a solution'.format(K)
            obj, clu = kclustering_fixed(points, K, restarts, TERROR, PAIRWISE, lord, j)
            objs[K] = obj
            clus[K] = clu
        log('Objective value for {}: {}'.format(K, objs[K]), level='INFO')

    MAXR = R
    compute(MAXR)
    
    while(R - L > 1):
        M = int(math.floor(float(R + L) / 2.0))
        assert M not in {L, R}, 'Median point is equal to boundaries but it cannot happen'
        compute(M)
        if objs[M] - objs[MAXR] > threshold:
            L = M
        else:
            R = M
            
    compute(L)
    compute(R)
    if L <= threshold:
        return clus[L]
    else:
        return clus[R]


def getord(lord):
    ## K-means minimizes SQUARE l2-norms while K-medians minimizes L1-norm
    if lord == 1: ## K-medians
        error = (lambda a, b : np.linalg.norm(a - b, ord=1))
        center = (lambda X : np.median(X, axis=0))
        pdist = (lambda X : scipy.spatial.distance.pdist(X, metric='cityblock'))
        cdist = (lambda X, Y : scipy.spatial.distance.cdist(X, Y, metric='cityblock'))
    elif lord == 2: ## K-means
        error = (lambda a, b : np.linalg.norm(a - b, ord=2)**2)
        center = (lambda X : np.mean(X, axis=0))
        pdist = (lambda X : scipy.spatial.distance.pdist(X, metric='sqeuclidean'))
        cdist = (lambda X, Y : scipy.spatial.distance.cdist(X, Y, metric='sqeuclidean'))
    else:
        raise ValueError('Order of l-norm distance must be either 1 or 2!')
    return error, center, pdist, cdist


def kclustering_fixed(points, K, restarts, TERROR, PAIRWISE, lord=1, j=1):
    with warnings.catch_warnings() as w:
        warnings.simplefilter("ignore")
        shared_points, shared_points_base = share_matrix(points)
        shared_pairwise, shared_pairwise_base = share_array(PAIRWISE)
        shared_clus, shared_clus_base = newshare_matrix(restarts, len(points))
    
    jobs = ((np.random.randint(low=0, high=2**10), x) for x, i in enumerate(range(restarts)))
    bar = ProgressBar(total=restarts, length=40, verbose=False)
    
    initargs = (points.shape[0], points.shape[1], K, lord, TERROR, shared_points, shared_pairwise, shared_clus)
    pool = Pool(processes=min(j, restarts), initializer=init_kclustering, initargs=initargs)
    progress = (lambda obj, it : bar.progress(advance=True, msg="Found clustering with obj: {} [Iterations: {}]".format(obj, it)))
    best = min(((obj, idx) for obj, idx, it in pool.imap_unordered(run_kclustering, jobs) if progress(obj, it)), key=(lambda x : x[0]))
    pool.close()
    pool.join()
    return best[0], shared_clus[best[1]]


def share_array(npdata):
    N = npdata.shape[0]
    shared_array_base = mp.Array(ctypes.c_double, N)
    shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
    shared_array[:] = npdata
    return shared_array, shared_array_base


def share_matrix(npdata):
    N, M = npdata.shape
    shared_matrix_base = mp.Array(ctypes.c_double, N * M)
    shared_matrix = np.ctypeslib.as_array(shared_matrix_base.get_obj())
    shared_matrix = shared_matrix.reshape(N, M)
    shared_matrix[:] = npdata
    return shared_matrix, shared_matrix_base


def newshare_matrix(N, M):
    shared_matrix_base = mp.Array(ctypes.c_double, N * M)
    shared_matrix = np.ctypeslib.as_array(shared_matrix_base.get_obj())
    shared_matrix = shared_matrix.reshape(N, M)
    return shared_matrix, shared_matrix_base


def init_kclustering(_N, _M, _K, _lord, _TERROR, _points, _pairwise, _clus):
    global N, M, K, error, center, pdist, cdist, TERROR, POINTS, PAIRWISE, CLUS
    N = _N
    M = _M
    K = _K
    error, center, pdist, cdist = getord(_lord)
    TERROR = _TERROR
    POINTS = _points
    PAIRWISE = _pairwise
    CLUS = _clus


def run_kclustering(job):
    seed, idx = job

    ## Utils
    np.random.seed(seed)
    randint = np.random.randint
    choice = np.random.choice
    lookup = (lambda i, j : PAIRWISE[indices_to_condensed(i, j, N)])

    ## Initialization
    centroids = []
    for i in range(K):
        if len(centroids) == 0:
            chosen = randint(N)
            probs = [lookup(x, chosen)**2 if x != chosen else 0.0 for x in xrange(N)]
        else:
            chosen = weighted_ichoice(probs)
            probs = [min(probs[x], lookup(x, chosen)**2) if x != chosen else 0.0 for x in xrange(N)]
        centroids.append(chosen)
    assert len(centroids) == K, 'Found less centroids {} than expected {}'.format(len(centroids), K)
    centroids = np.stack([POINTS[i] for i in centroids])

    ## Iterative process
    it = 0
    pre = None
    while pre is None or np.any(np.abs(centroids - pre) > 0.001):
        it += 1
        pre = centroids

        ## Assignment
        between  = cdist(POINTS, centroids)
        clu = [min((i for i in range(K)), key=(lambda i : between[x, i])) for x in xrange(N)]
        used = set(clu)

        ## Update centroids
        centerize = (lambda i : center(np.stack([p for x, p in enumerate(POINTS) if clu[x] == i])))
        centroids = np.stack([centerize(i) if i in used else np.zeros(M) for i in range(K)])

    CLUS[idx] = clu
        
    return sum(between[x, clu[x]] for x, p in enumerate(POINTS)) / TERROR, idx, it


def weighted_ichoice(weights):
    w = np.array(weights)
    wsum = np.sum(w, dtype=float)
    if wsum > 0:
        assert np.isfinite(wsum).all(), 'wsum distance contain NaN!\n{}'.format(wsum)
        assert np.isfinite(w).all(), 'w distance contain NaN!\n{}'.format(w)
        cs = np.cumsum(w) / wsum
        r = np.random.rand()
        return np.searchsorted(cs, r)
    else:
        return np.random.choice(np.arange(len(weights)), size=1)[0]

