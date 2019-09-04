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
    description = "CHISEL command to re-create the plots."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='?', default='calls/calls.tsv', type=str, help="Input file with inferred copy numbers (default: calls/calls.tsv)")
    parser.add_argument("-m", "--clonemap", required=False, type=str, default='clones/mapping.tsv', help="Clone map (default: not used, the cells will be clustered for plotting purposes)")
    parser.add_argument("-f", "--figformat", required=False, type=str, default='png', help="Format of output figures (default: png, the only other option is pdf)")
    parser.add_argument("-s", "--sample", required=False, type=int, default=20, help="Number of cells to sample (default: 20)")
    parser.add_argument("--excludenoisy", required=False, default=False, action='store_true', help="Exclude noisy cells from plots (default: False)")
    parser.add_argument("--gridsize", required=False, type=str, default='12,6', help="Grid dimenstions specified as comma-separated numbers (default: 12,6)")
    parser.add_argument("--plotsize", required=False, type=str, default='5,1.5', help="Plot dimenstions for RDR-BAF plots, specified as comma-separated numbers (default: 5,1.5)")
    parser.add_argument("--clussize", required=False, type=str, default='5,3', help="Grid dimenstions for clustered plots, specified as comma-separated numbers (default: 5,3)")
    parser.add_argument("--xmax", required=False, type=float, default=None, help="Maximum x-axis value (default: None)")
    parser.add_argument("--xmin", required=False, type=float, default=None, help="Minimum x-axis value (default: None)")
    parser.add_argument("--ymax", required=False, type=float, default=None, help="Maximum x-axis value (default: None)")
    parser.add_argument("--ymin", required=False, type=float, default=None, help="Minimum x-axis value (default: None)")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError('ERROR: input file {} does not exist!'.format(args.INPUT))
    if args.clonemap and not os.path.isfile(args.clonemap):
        raise ValueError('ERROR: the provided clone map does not exist!')
    if args.figformat not in ['pdf', 'png']:
        raise ValueError('ERROR: figure format must be either pdf or png!')
    if args.sample < 1:
        raise ValueError('ERROR: number of sampled cells must be positive!')

    return {
        'input' : args.INPUT,
        'clonemap' : args.clonemap,
        'format' : args.figformat,
        'sample' : args.sample,
        'nonoisy' : args.excludenoisy,
        'gridsize' : args.gridsize,
        'plotsize' : args.plotsize,
        'clussize' : args.clussize,
        'xmax' : args.xmax,
        'xmin' : args.xmin,
        'ymax' : args.ymax,
        'ymin' : args.ymin
    }


def main():
    log('Parsing and checking arguments', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories', level='PROGRESS')
    dplo = 'plots'
    if os.path.isdir(dplo):
        log("The plots sub-directory in the running directory already exists, results will be overwritten!", level='WARN')
    else:
        os.mkdir(dplo)

    log('Plotting', level='PROGRESS')
    os.chdir(dplo)
    up = (lambda f : os.path.join(os.pardir, f))
    cmd = 'python2.7 {} {} -m {} -f {} -s {}'
    cmd = cmd.format(os.path.join(src, 'Plotter.py'), up(args['input']), up(args['clonemap']), args['format'], args['sample'])
    cmd += ' --gridsize {}'.format(args['gridsize'])
    cmd += ' --plotsize {}'.format(args['plotsize'])
    cmd += ' --clussize {}'.format(args['clussize'])
    if args['nonoisy']:
        cmd += ' --excludenoisy'
    if args['xmax']:
        cmd += ' --xmax {}'.format(args['xmax'])
    if args['xmin']:
        cmd += ' --xmin {}'.format(args['xmin'])
    if args['ymax']:
        cmd += ' --ymax {}'.format(args['ymax'])
    if args['ymin']:
        cmd += ' --ymin {}'.format(args['ymin'])

    runcmd(cmd, './')
    os.chdir(os.pardir)


if __name__ == '__main__':
    main()
