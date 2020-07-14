#!/usr/bin/env python2.7

import sys, os
import argparse
import subprocess as sp
import multiprocessing as mp
import shlex
import datetime
import re

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'src')
if not os.path.isdir(src):
    raise ValueError("src directory not found in parent directory of bin i.e. {}, is anything been moved?".format(src))
sys.path.append(src)
from Utils import *


def parse_args():
    description = "CHISEL command to re-run the inference of allele- and haplotype-specific copy numbers, cell clustering, and plotting. This steps starts from estimated RDRs and BAFs."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='?', default='combo/combo.tsv', type=str, help="Input file with combined RDR and BAF per bin and per cell (default: combo/combo.tsv)")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-A","--sensitivity", required=False, type=float, default=1.0, help="Sensitivity of model selection for ploidy (default: 1, increase this parameter to lower sensitivity to noisy data, adjust this value (e.g. 2, 4, ..., 10, ...) to better deal with high-variance data (e.g. low coverage, small number of cells, low number of phased SNPs, etc...)")
    parser.add_argument("-P","--maxploidy", required=False, type=int, default=4, help="Maximum total copy number to consider for balanced cluster (default: 4, corresponding to a WGD)")
    parser.add_argument("-K","--upperk", required=False, type=int, default=100, help="Maximum number of bin clusters (default: 100, use 0 to consider maximum number of clusters)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError("Input file does not exist: {}".format(args.INPUT))
    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")
    if args.maxploidy < 3:
        raise ValueError("The maximum total copy number to consider for balanced cluster must be at least 2!")
    if args.upperk < 1:
        raise ValueError("The maximum number of clusters must be positive!")
    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")

    return {
        "INPUT" : args.INPUT,
        "rundir" : args.rundir,
        "sensitivity" : args.sensitivity,
        "maxploidy" : args.maxploidy,
        "upperk" : args.upperk,
        "seed" : args.seed,
        "jobs" : args.jobs
    }


def main():
    log('Parsing and checking arguments', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories', level='PROGRESS')
    dcal, dclo, dplo = setup(args, force=False)
    def get_comp(name):
        comp = os.path.join(src, name)
        if not os.path.isfile(comp):
            raise ValueError("{} not found in src directory of bin i.e. {}, is anything been moved?".format(name, src))
        return comp

    log('Calling', level='PROGRESS')
    cmd = 'python2.7 {} {} -A {} -P {} -K {} -j {}'
    cmd = cmd.format(get_comp('Caller.py'), args['INPUT'], args['sensitivity'], args['maxploidy'], args['upperk'], args['jobs'])
    if args['seed'] is not None:
        cmd += " --seed {}".format(args['seed'])
    runcmd(cmd, dcal, out='calls.tsv')
    calls = os.path.join(dcal, 'calls.tsv')

    log('Cloning', level='PROGRESS')
    cmd = 'python2.7 {} {}'
    cmd = cmd.format(get_comp('Cloner.py'), calls)
    if args['seed'] is not None:
        cmd += " --seed {}".format(args['seed'])
    runcmd(cmd, dclo, out='mapping.tsv')
    mapping = os.path.join(dclo, 'mapping.tsv')

    log('Plotting', level='PROGRESS')
    os.chdir(dplo)
    up = (lambda f : os.path.join(os.pardir, f))
    cmd = 'python2.7 {} {} -m {}'
    cmd = cmd.format(os.path.join(src, 'Plotter.py'), up(calls), up(mapping))
    runcmd(cmd, './')
    os.chdir(os.pardir)


def setup(args, force=True):
    dcal = os.path.join(args['rundir'], 'calls')
    if os.path.isdir(dcal):
        log("The calls sub-directory in the running directory already exists, results will be overwritten!", level='WARN')
    else:
        os.mkdir(dcal)

    dclo = os.path.join(args['rundir'], 'clones')
    if os.path.isdir(dclo):
        log("The clones sub-directory in the running directory already exists, results will be overwritten!", level='WARN')
    else:
        os.mkdir(dclo)

    dplo = os.path.join(args['rundir'], 'plots')
    if os.path.isdir(dplo):
        log("The plots sub-directory in the running directory already exists, results will be overwritten!\n", level='WARN')
    else:
        os.mkdir(dplo)

    return dcal, dclo, dplo


if __name__ == '__main__':
    main()
