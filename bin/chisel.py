#!/usr/bin/env python2.7

import sys, os
import argparse

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'src')
if not os.path.isdir(src):
    raise ValueError("src directory not found in parent directory of bin i.e. {}, is anything been moved?".format(src))
sys.path.append(src)
from Utils import *


def parse_args():
    description = "CHISEL command to run the complete pipeline starting from the 4 required data: (1) Barcoded single-cell BAM; (2) Matched-normal BAM; (3) Reference genome; (4) Phased VCF."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-t","--tumor", required=True, type=str, help="Barcoded single-cell BAM file")
    parser.add_argument("-n","--normal", required=True, type=str, help="Matched-normal BAM file")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-l","--listphased", type=str, required=True, help="Phased SNPs file (lines of heterozygous germline SNPs must contain either 0|1 or 1|0)")
    parser.add_argument("-b","--size", type=str, required=False, default="5Mb", help="Bin size, with or without \"kb\" or \"Mb\"")
    parser.add_argument("-k", "--blocksize", required=False, type=str, default="50kb", help="Size of the haplotype blocks (default: 50kb, use 0 to disable)")
    parser.add_argument("-c", "--chromosomes", type=str, required=False, default=' '.join(['chr{}'.format(i) for i in range(1, 23)]), help="Space-separeted list of chromosomes between apices (default: \"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\")")
    parser.add_argument("-m","--minreads", type=int, required=False, default=100000, help="Minimum number total reads to select cells (default: 100000)")
    parser.add_argument("-p","--maxploidy", required=False, type=int, default=4, help="Maximum total copy number to consider for balanced cluster (default: 4, corresponding to a WGD)")
    parser.add_argument("-K","--upperk", required=False, type=int, default=100, help="Maximum number of bin clusters (default: 100, use 0 to consider maximum number of clusters)")
    parser.add_argument("--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    args = parser.parse_args()

    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if not os.path.isfile(args.tumor):
        raise ValueError("Barcoded single-cell BAM file does not exist: {}".format(args.tumor))
    if not os.path.isfile(args.normal):
        raise ValueError("Matched-normal BAM file does not exist: {}".format(args.normal))
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome file does not exist: {}".format(args.reference))
    if not os.path.isfile(args.listphased):
        raise ValueError("Phased SNPs file does not exist: {}".format(args.listphased))
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")
    if args.minreads < 1:
        raise ValueError("The minimum number of reads must be positive!")
    if args.maxploidy < 3:
        raise ValueError("The maximum total copy number to consider for balanced cluster must be at least 2!")
    if args.upperk < 1:
        raise ValueError("The maximum number of clusters must be positive!")

    size = 0
    try:
        if args.size[-2:] == "kb":
            size = int(args.size[:-2]) * 1000
        elif args.size[-2:] == "Mb":
            size = int(args.size[:-2]) * 1000000
        else:
            size = int(args.size)
    except:
	raise ValueError("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!")

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
        "rundir" : args.rundir,
        "tumor" : args.tumor,
        "normal" : args.normal,
        "reference" : args.reference,
        "listphased" : args.listphased,
        "binsize" : size,
        "blocksize" : blocksize,
        "chromosomes" : args.chromosomes,
        "minreads" : args.minreads,
        "bcftools" : args.bcftools,
        "samtools" : args.samtools,
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
    dbaf, drdr, dcom, dcal, dclo, dplo = setup(args)
    def get_comp(name):
        comp = os.path.join(src, name)
        if not os.path.isfile(comp):
            raise ValueError("{} not found in src directory of bin i.e. {}, is anything been moved?".format(name, src))
        return comp

    log('Computing RDRs', level='PROGRESS')
    cmd = 'python2.7 {} -n {} -t {} -r {} -b {} -m {} -j {} --outdir {}'
    cmd = cmd.format(get_comp('RDREstimator.py'), args['normal'], args['tumor'], args['reference'], args['binsize'], args['minreads'], args['jobs'], drdr)
    if args['samtools'] is not None:
        cmd += " -s {}".format(args['samtools'])
    runcmd(cmd, drdr, out='rdr.tsv')
    lcel = os.path.join(drdr, 'total.tsv')
    rdr = os.path.join(drdr, 'rdr.tsv')

    log('Computing BAFs', level='PROGRESS')
    cmd = 'python2.7 {} -n {} -t {} -r {} -j {} -c {} -l {}'
    cmd = cmd.format(get_comp('BAFEstimator.py'), args['normal'], args['tumor'], args['reference'], args['jobs'], lcel, args['listphased'])
    if args['samtools'] is not None:
        cmd += " -s {}".format(args['samtools'])
    if args['bcftools'] is not None:
        cmd += " -b {}".format(args['bcftools'])
    runcmd(cmd, dbaf, out='baf.tsv')
    baf = os.path.join(dbaf, 'baf.tsv')

    log('Combining RDRs and BAFs', level='PROGRESS')
    cmd = 'python2.7 {} -r {} -b {} -j {} -k {} -l {}'
    cmd = cmd.format(get_comp('Combiner.py'), rdr, baf, args['jobs'], args['blocksize'], lcel)
    if args['seed'] is not None:
        cmd += " --seed {}".format(args['seed'])
    runcmd(cmd, dcom, out='combo.tsv')
    com = os.path.join(dcom, 'combo.tsv')

    log('Calling', level='PROGRESS')
    cmd = 'python2.7 {} {} -p {} -K {} -j {}'
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
    cmd = cmd.format(os.path.join(src, 'Plotter.py'), up(args['input']), up(mapping))
    runcmd(cmd, './')
    os.chdir(os.pardir)


def setup(args):
    dbaf = os.path.join(args['rundir'], 'baf')
    if os.path.isdir(dbaf):
        raise ValueError("The baf directory in the running directory already exists! Please remove it or change running directory to avoid overwriting!")
    else:
        os.mkdir(dbaf)

    drdr = os.path.join(args['rundir'], 'rdr')
    if os.path.isdir(drdr):
        raise ValueError("The rdr directory in the running directory already exists! Please remove it or change running directory to avoid overwriting!")
    else:
        os.mkdir(drdr)

    dcom = os.path.join(args['rundir'], 'combo')
    if os.path.isdir(dcom):
        raise ValueError("The combo directory in the running directory already exists! Please remove it or change running directory to avoid overwriting!")
    else:
        os.mkdir(dcom)

    dcal = os.path.join(args['rundir'], 'calls')
    if os.path.isdir(dcal):
        raise ValueError("The calls directory in the running directory already exists! Please remove it or change running directory to avoid overwriting!")
    else:
        os.mkdir(dcal)

    dclo = os.path.join(args['rundir'], 'clones')
    if os.path.isdir(dclo):
        raise ValueError("The clones directory in the running directory already exists! Please remove it or change running directory to avoid overwriting!")
    else:
        os.mkdir(dclo)

    dplo = os.path.join(args['rundir'], 'plots')
    if os.path.isdir(dplo):
        raise ValueError("The plots directory in the running directory already exists! Please remove it or change running directory to avoid overwriting!")
    else:
        os.mkdir(dplo)

    return dbaf, drdr, dcom, dcal, dclo, dplo


if __name__ == '__main__':
    main()
