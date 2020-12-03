#!/usr/bin/env python2.7

import os
import argparse
import subprocess as sp
import multiprocessing as mp
import shlex
import datetime
import re

import chisel

src = os.path.dirname(chisel.__file__)
from ..Utils import *


def parse_args():
    description = "CHISEL command to run the pipeline starting from inferred copy numbers."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='?', default='calls/calls.tsv', type=str, help="Input file with combined RDR and BAF per bin and per cell")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-f", "--maxdiff", required=False, type=float, default=0.06, help="Maximum haplotype-specific distance between the genome of cells in the same clone (default: 0.06, when -1 is chosen the maximum cluster method of SciPy is used)")
    parser.add_argument("-s", "--minsize", required=False, type=int, default=14, help="Minimum number of cells in a subpopulation to define a clone (default: 14)")
    parser.add_argument("-K","--upperk", required=False, type=int, default=100, help="Maximum number of bin clusters (default: 100, use 0 to consider maximum number of clusters)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError("Input file does not exist: {}".format(args.INPUT))
    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")
    if (args.maxdiff < 0.0 and args.maxdiff != 1.0) or args.maxdiff > 1.0:
        raise ValueError("Maximum distance must be in [0, 1] or equal to -1!")
    if args.minsize < 0:
        raise ValueError("Minimum number of cells in a clone must be positive!")
    if args.upperk < 1:
        raise ValueError("The maximum number of clusters must be positive!")

    return {
        "INPUT" : args.INPUT,
        "rundir" : args.rundir,
        "maxdiff" : args.maxdiff,
        "minsize" : args.minsize,
        "upperk" : args.upperk,
        "seed" : args.seed
    }


def main():
    log('Parsing and checking arguments', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories', level='PROGRESS')
    dclo, dplo = setup(args, force=False)
    def get_comp(name):
        comp = os.path.join(src, name)
        if not os.path.isfile(comp):
            raise ValueError("{} not found in src directory of bin i.e. {}, is anything been moved?".format(name, src))
        return comp

    log('Cloning', level='PROGRESS')
    cmd = 'python2.7 {} {} -f {} -s {}'
    cmd = cmd.format(get_comp('Cloner.py'), args['INPUT'], args['maxdiff'], args['minsize'])
    if args['seed'] is not None:
        cmd += " --seed {}".format(args['seed'])
    runcmd(cmd, dclo, out='mapping.tsv')
    mapping = os.path.join(dclo, 'mapping.tsv')

    log('Plotting', level='PROGRESS')
    os.chdir(dplo)
    up = (lambda f : os.path.join(os.pardir, f))
    cmd = 'python2.7 {} {} -m {}'
    cmd = cmd.format(os.path.join(src, 'Plotter.py'), up(args['INPUT']), up(mapping))
    runcmd(cmd, './')
    os.chdir(os.pardir)


def setup(args, force=True):
    dclo = os.path.join(args['rundir'], 'clones')
    if os.path.isdir(dclo):
        log("The clones sub-directory in the running directory already exists, results will be overwritten!", level='WARN')
    else:
        os.mkdir(dclo)

    dplo = os.path.join(args['rundir'], 'plots')
    if os.path.isdir(dplo):
        log("The plots sub-directory in the running directory already exists, results will be overwritten!", level='WARN')
    else:
        os.mkdir(dplo)

    return dclo, dplo


if __name__ == '__main__':
    main()
