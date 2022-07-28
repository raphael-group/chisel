#!/usr/bin/env python2.7

import os, sys
os.environ["OMP_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
import argparse
import shlex
import shutil
import glob
import re
import traceback
import subprocess as sp
import multiprocessing as mp

import numpy as np

import chisel
src = os.path.dirname(chisel.__file__)
from ..Utils import *
from chisel.RDREstimator import get_bins


def parse_args():
    description = "CHISEL command to only compute RDRs per cell."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-t","--tumor", required=True, type=str, help="Barcoded single-cell BAM file")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-b","--size", type=str, required=False, default="250kb", help="Bin size, with or without \"kb\" or \"Mb\" (default: 250kb)")
    parser.add_argument("-c", "--chromosomes", type=str, required=False, default=' '.join(['chr{}'.format(i) for i in range(1, 23)]), help="Space-separeted list of chromosomes between apices (default: \"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\")")
    parser.add_argument("-m","--minreads", type=int, required=False, default=100000, help="Minimum number total reads to select cells (default: 100000)")
    parser.add_argument("--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("--art", required=False, default=None, type=str, help="Path to the directory to \"art_illumina\" executable (default: in $PATH)")
    parser.add_argument("--bwa", required=False, default=None, type=str, help="Path to the directory to \"bwa\" executable (default: in $PATH)")
    parser.add_argument("--cellprefix", type=str, required=False, default='CB:Z:', help="Prefix of cell barcode field in SAM format (default: CB:Z:)")
    parser.add_argument("--cellsuffix", type=str, required=False, default=None, help="Suffix of cell barcode field in SAM format (default: none)")
    parser.add_argument("--simcov", required=False, type=float, default=.2, help="Sequencing fold coverage of simulated normal BAM file (default: 0.2)")
    parser.add_argument("--binstats", required=False, type=int, default=None, help="Number of bins to sample per chromosome to estimate sequencing stats (default: all are used, fix a number for improving speed)")
    parser.add_argument("--keeptmpdir", required=False, default=False, action='store_true', help="Do not erase temporary directory (default: False)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if not os.path.isfile(args.tumor):
        raise ValueError("Barcoded single-cell BAM file does not exist: {}".format(args.tumor))
    if args.seed and args.seed < 1:
        raise ValueError("The random seed  must be positive!")
    if args.minreads < 1:
        raise ValueError("The minimum number of reads must be positive!")
    if args.simcov <= 0.0:
        raise ValueError("The sequencing coverage of simulated normal must be >= 0.0!")
    if args.binstats is not None and args.binstats <= 0:
        raise ValueError("The number of bins for sequencing stats must be >= 0.0!")

    if not os.path.isfile(args.reference):
        raise ValueError(error("Reference genome file does not exist: {}".format(args.reference)))
    refidx = ['{}.{}'.format(args.reference, ix) for ix in ['amb', 'ann', 'bwt', 'pac', 'sa']]
    if not all(os.path.isfile(f) for f in refidx):
        raise ValueError(error("Some of the BWA index files are missing, please make sure these are available and generated through the command \n\t``bwa index {}''.\n Expected files are: {}".format(args.reference, '\n'.join(refidx))))

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

    art = args.art
    if not art:
        art = "art_illumina"
    if which(art) is None:
        raise ValueError(error("art_illumina has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} bwa\n\nOtherwise, please provide with the flag `--art` the full path to the directory containing art_illumina exacutable."))

    bwa = args.bwa
    if not bwa:
        bwa = "bwa"
    if which(bwa) is None and not args.barcodeonly:
        raise ValueError(error("bwa has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} bwa\n\nOtherwise, please provide with the flag `--bwa` the full path to the directory containing bwa exacutable."))


    return {
        "rundir" : os.path.abspath(args.rundir),
        "tumor" : os.path.abspath(args.tumor),
        "reference" : os.path.abspath(args.reference),
        "binsize" : size,
        "chromosomes" : args.chromosomes,
        "minreads" : args.minreads,
        "bcftools" : bcftools,
        "samtools" : samtools,
        "bwa" : bwa,
        "art" : art,
        "cellprefix" : args.cellprefix,
        "cellsuffix" : args.cellsuffix,
        "simcov" : args.simcov,
        "binstats" : args.binstats,
        "keeptmpdir" : args.keeptmpdir,
        "seed" : args.seed,
        "jobs" : args.jobs
    }


def main():
    log('Parsing and checking arguments', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')

    log('Setting directories', level='PROGRESS')
    ddata, dbaf, drdr, dcom, dcal, dclo, dplo = setup(args)
    def get_comp(name):
        comp = os.path.join(src, name)
        if not os.path.isfile(comp):
            raise ValueError("{} not found in src directory of bin i.e. {}, is anything been moved?".format(name, src))
        return comp

    log('Simulating DNA sequencing reads for mappability correction', level='PROGRESS')
    args['normal'] = os.path.join(ddata, 'simulatednormal.bam')
    simulate_sequencing(args)

    log('Computing RDRs', level='PROGRESS')
    cmd = 'python2.7 {} -n {} -t {} -r {} -b {} -m {} -j {} -c \"{}\" --outdir {}'
    cmd = cmd.format(get_comp('RDREstimator.py'), args['normal'], args['tumor'], args['reference'], args['binsize'], args['minreads'], args['jobs'], args['chromosomes'], drdr)
    if args['samtools'] is not None:
        cmd += " -s {}".format(args['samtools'])
    cmd += " --cellprefix {}".format(args['cellprefix'])
    if args['cellsuffix'] is not None:
        cmd += " --cellsuffix {}".format(args['cellsuffix'])
    runcmd(cmd, drdr, out='rdr.tsv')
    lcel = os.path.join(drdr, 'total.tsv')
    rdr = os.path.join(drdr, 'rdr.tsv')


def setup(args):
    if any(os.path.isdir(os.path.join(args['rundir'], x)) for x in ['data', 'baf', 'rdr', 'combo', 'calls', 'clones', 'plots']):
        log('Some of the working folders already exist in the running directory and content will be overwritten, please interrupt the process if this was not intended.', level='WARN')

    ddata = os.path.join(args['rundir'], 'data')
    if not os.path.isdir(ddata):
        os.mkdir(ddata)
    
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

    return ddata, dbaf, drdr, dcom, dcal, dclo, dplo


def simulate_sequencing(args):
    log('Retrieving sequencing stats', level='INFO')
    bins = get_bins(args['reference'], args['chromosomes'].split(), args['binsize'], bams=[args['tumor']], samtools=args['samtools'])
    if args['binstats'] is not None:
        bins = {c : [bins[c][x] for x in np.random.choice(len(bins[c]), min(args['binstats'], len(bins[c])))] for c in bins}
    jobs = [(chro, cbin, args) for chro in bins for cbin in bins[chro]]
    bar = ProgressBar(total=len(jobs), length=30, verbose=False)
    pool = mp.Pool(processes=min(args['jobs'], len(jobs)))
    progress = (lambda e : bar.progress(advance=True, msg="Stats retrieved for {}:{}-{}".format(*e)))
    bamstats = [m for e, m in pool.imap_unordered(get_bam_stats, jobs) if progress(e)]
    bamstats = list(filter(lambda v : v['read'] > 0.0 or v['fragment'] > 0.0, bamstats))
    
    avg = (lambda L : sum(L) / float(len(L)))
    if any(s['fragment'] > 0 for s in bamstats):
        bamstats = [s for s in bamstats if s['read'] > 0.0 or s['fragment'] > 0.0]
        read = int(round(avg([s['read'] for s in bamstats])))
        fragment = int(round(avg([s['fragment'] for s in bamstats])))
        sdfragment = int(round(avg([s['sd'] for s in bamstats])))
        log('BAM has been identified as paired-end sequencing with read length {} and fragment size {} (sd: {})'.format(read, fragment, sdfragment), level='INFO')
        sequencing_paired(args, read, fragment, sdfragment)
    else:
        bamstats = [s for s in bamstats if s['read'] > 0.0]
        read = int(round(avg([s['read'] for s in bamstats])))
        log('BAM has been identified as paired-end sequencing with read length {}'.format(read), level='INFO')
        sequencing_single(args, read)
    assert os.path.isfile(args['normal'])
    log('Simulated normal sequencing reads are written at: {}'.format(args['normal']), level='INFO')


def get_bam_stats(dt):
    chro, bins, args = dt
    cmd_sam = '{} stats {} {}:{}-{}'.format(args['samtools'], args['tumor'], chro, bins[0], bins[1])
    cmd_gre = 'grep -e "average length" -e "insert size average" -e "insert size standard deviation"'
    sam = sp.Popen(shlex.split(cmd_sam), stdout=sp.PIPE, stderr=sp.PIPE)
    gre = sp.Popen(shlex.split(cmd_gre), stdin=sam.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = gre.communicate()
    if gre.returncode != 0:
        raise ValueError(error('Merging failed with messages:\n{}\n{}\n'.format(stdout, stderr)))
    result = stdout.strip().split('\n')
    if len(result) != 3:
        raise ValueError(error('More or less than expected BAM statistics have been found:\n{}\n'.format('\n'.join(result))))
    process = (lambda r : 'read' if 'average length' in r else ('fragment' if 'insert size average' in r else 'sd'))
    check = (lambda v : float(v) if v.replace('.', '', 1).isdigit() else 0)
    return (chro, bins[0], bins[1]), {process(rec) : check(rec.split(':')[-1].strip()) for rec in result}


def sequencing_paired(args, read, fragment, sdfragment):
    tmpdir = os.path.join(os.path.dirname(args['normal']), '_TMPDIR')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    log('Simulating sequencing reads', level='INFO')
    pre_fastq = os.path.join(tmpdir, 'normal')
    profile = 'HS20' if abs(100 - read) < abs(125 - read) else 'HS25'
    cmd = '{} -ss {} -na -i {} -p -l {} -f {} -m {} -s {} -o {}'.format(args['art'], profile, args['reference'], read, args['simcov'], fragment, sdfragment, pre_fastq)
    with open(os.path.join(tmpdir, 'art_illumina.log'), 'w') as simlog:
        proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=simlog)
        stdout, stderr = proc.communicate()
    fastqs = ('{}1.fq'.format(pre_fastq), '{}2.fq'.format(pre_fastq))
    if proc.returncode != 0 or not (os.path.exists(fastqs[0]) and os.path.exists(fastqs[1])):
        raise ValueError(error('ART Illumina simuation of sequencing reads failed:\n{}\n{}\n'.format(stdout, stderr)))
    
    log('Aligning simulated sequencing reads', level='INFO')
    jobs = max(1, int(round(args['jobs'] / 2.0)))
    cmd_bwa = '{} mem -M -t {} {} {} {}'.format(args['bwa'], jobs, args['reference'], fastqs[0], fastqs[1])
    cmd_sor = '{} sort - -Obam -o {} -T {} -@ {}'.format(args['samtools'], args['normal'], tmpdir, max(1, args['jobs'] - jobs))
    bwa = sp.Popen(shlex.split(cmd_bwa), stdout=sp.PIPE, stderr=sp.PIPE)
    sor = sp.Popen(shlex.split(cmd_sor), stdin=bwa.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sor.communicate()
    if sor.returncode != 0 or not os.path.exists(args['normal']):
        raise ValueError(error('Alignment failed with messages:\n{}\n{}\n'.format(stdout, stderr)))
    
    log('Indexing simulated sequencing reads', level='INFO')
    cmd_sam = '{} index {}'.format(args['samtools'], args['normal'])
    sam = sp.Popen(shlex.split(cmd_sam), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sam.communicate()
    fastqs = ('{}1.fq'.format(pre_fastq), '{}2.fq'.format(pre_fastq))
    if sor.returncode != 0:
        raise ValueError(error('ART Illumina simuation of sequencing reads failed:\n{}\n{}\n'.format(stdout, stderr)))
    
    if not args['keeptmpdir']:
        shutil.rmtree(tmpdir)
    return 


def sequencing_single(args, read):
    tmpdir = os.path.join(os.path.dirname(args['normal']), '_TMPDIR')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    log('Simulating sequencing reads', level='INFO')
    pre_fastq = os.path.join(tmpdir, 'normal')
    profile = 'HS20' if abs(100 - read) < abs(125 - read) else 'HS25'
    cmd = '{} -ss {} -na -i {} -l {} -f {} -o {}'.format(args['art'], profile, args['reference'], read, args['simcov'], pre_fastq)
    with open(os.path.join(tmpdir, 'art_illumina.log'), 'w') as simlog:
        proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=simlog)
        stdout, stderr = proc.communicate()
    fastq = '{}.fq'.format(pre_fastq)
    if proc.returncode != 0 or not os.path.exists(fastq):
        raise ValueError(error('ART Illumina simuation of sequencing reads failed:\n{}\n{}\n'.format(stdout, stderr)))
    
    log('Aligning simulated sequencing reads', level='INFO')
    jobs = max(1, int(round(args['jobs'] / 2.0)))
    cmd_bwa = '{} mem -M -t {} {} {}'.format(args['bwa'], jobs, args['reference'], fastq)
    cmd_sor = '{} sort - -Obam -o {} -T {} -@ {}'.format(args['samtools'], args['normal'], tmpdir, max(1, args['jobs'] - jobs))
    bwa = sp.Popen(shlex.split(cmd_bwa), stdout=sp.PIPE, stderr=sp.PIPE)
    sor = sp.Popen(shlex.split(cmd_sor), stdin=bwa.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sor.communicate()
    if sor.returncode != 0 or not os.path.exists(args['normal']):
        raise ValueError(error('Alignment failed with messages:\n{}\n{}\n'.format(stdout, stderr)))
    
    log('Indexing simulated sequencing reads', level='INFO')
    cmd_sam = '{} index {}'.format(args['samtools'], args['normal'])
    sam = sp.Popen(shlex.split(cmd_sam), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sam.communicate()
    fastqs = ('{}1.fq'.format(pre_fastq), '{}2.fq'.format(pre_fastq))
    if sor.returncode != 0:
        raise ValueError(error('ART Illumina simuation of sequencing reads failed:\n{}\n{}\n'.format(stdout, stderr)))
    
    if not args['keeptmpdir']:
        shutil.rmtree(tmpdir)
    return 


if __name__ == '__main__':
    main()
