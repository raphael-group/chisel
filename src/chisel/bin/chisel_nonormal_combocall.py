#!/usr/bin/env python2.7

import os, sys
os.environ["OMP_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
import argparse
from subprocess import Popen
import chisel

src = os.path.dirname(chisel.__file__)
from ..Utils import *


def parse_args():
    description = "CHISEL command to run the complete pipeline starting from RDRs and BAFs for one or multiple samples from previously executions of CHISEL or CHISEl preprocess."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", type=str, nargs='+', help="One or multiple CHISEL directory runs for different samples from which to combine RDRs and BAFs")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("--names", required=False, default=None, type=str, nargs='+', help="Sample names when combining multiple samples (default: idx used)")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-k", "--blocksize", required=False, type=str, default="50kb", help="Size of the haplotype blocks (default: 50kb, use 0 to disable)")
    parser.add_argument("-c", "--chromosomes", type=str, required=False, default=' '.join(['chr{}'.format(i) for i in range(1, 23)]), help="Space-separeted list of chromosomes between apices (default: \"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\")")
    parser.add_argument("-p","--maxploidy", required=False, type=int, default=4, help="Maximum total copy number to consider for balanced cluster (default: 4, corresponding to a WGD)")
    parser.add_argument("-K","--upperk", required=False, type=int, default=100, help="Maximum number of bin clusters (default: 100, use 0 to consider maximum number of clusters)")
    parser.add_argument("--minimumsnps", required=False, type=float, default=0.08, help="Minimum SNP density per kb (default: 0.08)")
    parser.add_argument("--missingsnps", required=False, type=str, default="10,0", help="A,B counts for genomic bins without minimum minimum SNP density (default: 10,0 i.e. BAF=0)")
    parser.add_argument("--nogccorr", required=False, default=False, action='store_true', help="Disable correction for GC bias (default: enabled)")
    parser.add_argument("--nophasecorr", required=False, default=False, action='store_true', help="Disable correction for given phasing bias (default: enabled)")
    parser.add_argument("--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("--cellprefix", type=str, required=False, default='CB:Z:', help="Prefix of cell barcode field in SAM format (default: CB:Z:)")
    parser.add_argument("--cellsuffix", type=str, required=False, default=None, help="Suffix of cell barcode field in SAM format (default: none)")
    parser.add_argument("--simcov", required=False, type=float, default=2, help="Sequencing fold coverage of simulated normal BAM file (default: 2)")
    parser.add_argument("--binstats", required=False, type=int, default=None, help="Number of bins to sample per chromosome to estimate sequencing stats (default: all are used, fix a number for improving speed)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    args = parser.parse_args()

    for indir in args.INPUT:
        if not os.path.isdir(indir):
            raise ValueError("Input directory does not exists: {}".format(indir))
        rdr_file = os.path.join(indir, 'rdr', 'rdr.tsv')
        if not os.path.isfile(rdr_file):
            raise ValueError("Input directory does not contain RDR file: {}".format(rdr_file))
        tot_file = os.path.join(indir, 'rdr', 'total.tsv')
        if not os.path.isfile(tot_file):
            raise ValueError("Input directory does not contain Total read file: {}".format(tot_file))
        baf_file = os.path.join(indir, 'baf', 'baf.tsv')
        if not os.path.isfile(baf_file):
            raise ValueError("Input directory does not contain BAF file: {}".format(baf_file))
    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")
    if args.maxploidy < 3:
        raise ValueError("The maximum total copy number to consider for balanced cluster must be at least 2!")
    if args.upperk < 1:
        raise ValueError("The maximum number of clusters must be positive!")
    if args.minimumsnps < 0.0:
        raise ValueError("The minimum SNP density must be >= 0.0!")
    if args.simcov <= 0.0:
        raise ValueError("The sequencing coverage of simulated normal must be >= 0.0!")
    if args.binstats is not None and args.binstats <= 0:
        raise ValueError("The number of bins for sequencing stats must be >= 0.0!")

    if not os.path.isfile(args.reference):
        raise ValueError(error("Reference genome file does not exist: {}".format(args.reference)))
    refidx = ['{}.{}'.format(args.reference, ix) for ix in ['amb', 'ann', 'bwt', 'pac', 'sa']]
    if not all(os.path.isfile(f) for f in refidx):
        raise ValueError(error("Some of the BWA index files are missing, please make sure these are available and generated through the command \n\t``bwa index {}''.\n Expected files are: {}".format(args.reference, '\n'.join(refidx))))

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

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")

    bcftools = args.bcftools
    if not bcftools:
        bcftools = "bcftools"
    if which(bcftools) is None:
        raise ValueError("bcftools has not been found or is not executable!")

    samtools = args.samtools
    if not samtools:
        samtools = "samtools"
    if which(samtools) is None:
        raise ValueError("samtools has not been found or is not executable!")

    return {
        "input" : args.INPUT,
        "reference" : os.path.abspath(args.reference),
        "names" : args.names,
        "rundir" : args.rundir,
        "blocksize" : blocksize,
        "chromosomes" : args.chromosomes,
        "phasecorr" : not args.nophasecorr,
        "bcftools" : bcftools,
        "samtools" : samtools,
        "maxploidy" : args.maxploidy,
        "upperk" : args.upperk,
        'minimumsnps' : args.minimumsnps,
        'missingsnps' : args.missingsnps,
        "cellprefix" : args.cellprefix,
        "cellsuffix" : args.cellsuffix,
        "gccorr" : not args.nogccorr,
        "phasecorr" : not args.nophasecorr,
        "simcov" : args.simcov,
        "binstats" : args.binstats,
        "seed" : args.seed,
        "jobs" : args.jobs
    }


def main():
    log('Parsing and checking arguments', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories', level='PROGRESS')
    dbaf, drdr, dcom, dcal, dclo, dplo = setup(args)
    def get_comp(name):
        comp = os.path.join(src, name)
        if not os.path.isfile(comp):
            raise ValueError("{} not found in src directory of bin i.e. {}, is anything been moved?".format(name, src))
        return comp

    lcel = os.path.join(drdr, 'total.tsv')
    if os.path.isfile(lcel):
        raise ValueError("Total read file {} already exists, please remove it or it'd get overwritten!".format(lcel))
    rdr = os.path.join(drdr, 'rdr.tsv')
    if os.path.isfile(rdr):
        raise ValueError("RDR file {} already exists, please remove it or it'd get overwritten!".format(rdr))
    baf = os.path.join(dbaf, 'baf.tsv')
    if os.path.isfile(baf):
        raise ValueError("BAF file {} already exists, please remove it or it'd get overwritten!".format(baf))
    
    log('Aggregating previously-computed RDRs and BAFs', level='PROGRESS')
    aggregate(rdr, lcel, baf, args['input'], args['names'])

    log('Combining RDRs and BAFs', level='PROGRESS')
    cmd = 'python2.7 {} -r {} -b {} -j {} -k {} -l {} --minimumsnps {} --missingsnps {}'
    cmd = cmd.format(get_comp('Combiner.py'), rdr, baf, args['jobs'], args['blocksize'], lcel, args['minimumsnps'], args['missingsnps'])
    if args['gccorr']:
        cmd += " --gccorr {}".format(args['reference'])
    if not args['phasecorr']:
        cmd += " --nophasecorr"
    if args['seed'] is not None:
        cmd += " --seed {}".format(args['seed'])
    runcmd(cmd, dcom, out='combo.tsv')
    com = os.path.join(dcom, 'combo.tsv')

    log('Calling', level='PROGRESS')
    cmd = 'python2.7 {} {} -P {} -K {} -j {}'
    cmd = cmd.format(get_comp('Caller.py'), com, args['maxploidy'], args['upperk'], args['jobs'])
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


def aggregate(rdr, lcel, baf, input_dirs, names):
    if names is None or len(names) != len(input_dirs):
        names = list(range(len(input_dirs)))
    for indir, name in zip(input_dirs, names):
        log('Aggregating RDRs for {} with name {}'.format(indir, name), level='INFO')
        Popen("awk -v OFS='\\t' '{}' {} >> {}".format('{}print $1,$2,$3,"{}_"$4,$5,$6,$7{}'.format('{', name, '}'),
                                                     os.path.join(indir, 'rdr', 'rdr.tsv'),
                                                     rdr), shell=True).communicate() 
    for indir, name in zip(input_dirs, names):
        log('Aggregating Total reads for {} with name {}'.format(indir, name), level='INFO')
        Popen("awk -v OFS='\\t' '{}' {} >> {}".format('{}print "{}_"$1,$2{}'.format('{', name, '}'),
                                                     os.path.join(indir, 'rdr', 'total.tsv'),
                                                     lcel), shell=True).communicate() 
    for indir, name in zip(input_dirs, names):
        log('Aggregating BAFs for {} with name {}'.format(indir, name), level='INFO')
        Popen("awk -v OFS='\\t' '{}' {} >> {}".format('{}print $1,$2,"{}_"$3,$4,$5{}'.format('{', name, '}'),
                                                     os.path.join(indir, 'baf', 'baf.tsv'),
                                                     baf), shell=True).communicate() 
    return


def setup(args):
    if any(os.path.isdir(os.path.join(args['rundir'], x)) for x in ['baf', 'rdr', 'combo', 'calls', 'clones', 'plots']):
        log('Some of the working folders already exist in the running directory and content will be overwritten, please interrupt the process if this was not intended.', level='WARN')

    dbaf = os.path.join(args['rundir'], 'baf')
    if not os.path.isdir(dbaf):
        os.mkdir(dbaf)

    drdr = os.path.join(args['rundir'], 'rdr')
    if not os.path.isdir(drdr):
        os.mkdir(drdr)

    dcom = os.path.join(args['rundir'], 'combo')
    if not os.path.isdir(dcom):
        os.mkdir(dcom)

    dcal = os.path.join(args['rundir'], 'calls')
    if not os.path.isdir(dcal):
        os.mkdir(dcal)

    dclo = os.path.join(args['rundir'], 'clones')
    if not os.path.isdir(dclo):
        os.mkdir(dclo)

    dplo = os.path.join(args['rundir'], 'plots')
    if not os.path.isdir(dplo):
        os.mkdir(dplo)

    return dbaf, drdr, dcom, dcal, dclo, dplo


if __name__ == '__main__':
    main()
