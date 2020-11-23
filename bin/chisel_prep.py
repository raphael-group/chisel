#!/usr/bin/env python2.7

import sys, os
import argparse
import shlex
import shutil
import glob
import subprocess as sp
import multiprocessing as mp

from itertools import combinations_with_replacement
from heapq import nlargest

import numpy as np

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'src')
if not os.path.isdir(src):
    raise ValueError("src directory not found in parent directory of bin i.e. {}, is anything been moved?".format(src))
sys.path.append(src)
from Utils import *


def parse_args():
    description = "CHISEL command to create a barcoded BAM file from single-cell FASTQs, single-cell BAMs, or a `RG:Z:`-barcoded BAM files without `CB:Z:` field."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='+', type=str, help="Input FASTQs, BAMs, or a single BAM file with different behaviors: ................................. (1) FASTQs -- specified in a directory DIR as `DIR/*.fastq` or `DIR/*.fastq.gz` -- will be barcoded and aligned into a barcoded BAM file; .............. (2) BAMs -- specified in a directory DIR as `DIR/*.bam` -- will be barcoded and aligned into a barcoded BAM file; ................................. (3) a single BAM file with barcodes in the field `RG:Z:` will be converted into a barcoded BAM file with the additional `CB:Z:` field.")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-o","--output", required=False, default='barcodedcells.bam', type=str, help="Output name in running directory (default: barcodedcells.bam)")    
    parser.add_argument("--noduplicates", required=False, default=False, action='store_true', help="Do not perform marking duplicates and recalibration with Picard tools (default: False)")
    parser.add_argument("--keeptmpdir", required=False, default=False, action='store_true', help="Do not erase temporary directory (default: False)")
    parser.add_argument("--barcodelength", required=False, type=int, default=12, help="Length of barcodes (default: 12)")
    parser.add_argument("--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable (default: in $PATH)")
    parser.add_argument("--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable (default: in $PATH)")
    parser.add_argument("--bwa", required=False, default=None, type=str, help="Path to the directory to \"bwa\" executable (default: in $PATH)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    args = parser.parse_args()
    
    if args.seed is not None:
        np.random.seed(args.seed)

    def ispaired(L):
        D = defaultdict(lambda : [])
        map(lambda l : D[l[:-1]].append(l[-1]), L)
        return all(D[d] == [1, 2] or D[d] == [1, 2] for d in D)
    
    inputs = args.INPUT.split()
    if all(f[:-6] == '.fastq' for f in inputs):
        mode = 'q2' if ispaired(map(lambda f : f[:-6], inputs)) else 'q1'
    elif all(f[:-9] == '.fastq.gz' for f in inputs):
        mode = 'z2' if ispaired(map(lambda f : f[:-9], inputs)) else 'z1'
    elif all(f[:-4] == '.bam' for f in inputs):
        mode = 'B' if len(inputs) == 1 else 'b'
    else:
        raise ValueError("Input files are of wrong format or mixed formats")    

    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome file does not exist: {}".format(args.reference))

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError("The number of jobs must be positive!")

    bcftools = args.bcftools
    if not bcftools:
        bcftools = "bcftools"
    if which(bcftools) is None:
        raise ValueError("bcftools has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} bcftools\n\nOtherwise, please provide with the flag `--bcftools` the full path to the directory containing bcftools exacutable.")

    samtools = args.samtools
    if not samtools:
        samtools = "samtools"
    if which(samtools) is None:
        raise ValueError("samtools has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} samtools\n\nOtherwise, please provide with the flag `--samtools` the full path to the directory containing samtools exacutable.")

    bwa = args.bwa
    if not bwa:
        bwa = "bwa"
    if which(bwa) is None:
        raise ValueError("bwa has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} bwa\n\nOtherwise, please provide with the flag `--bwa` the full path to the directory containing bwa exacutable.")

    if mode in ['z1', 'z2']:
        if which('gzip') is None:
            raise ValueError("gzip has not been found or is not executable but is required for provided FASTQ files!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -n ${ENV} gzip\n\nOtherwise, please make it available in PATH.")
    
    output = os.path.basename(args.output if args.output[-4] == '.bam' else '{}.bam'.format(args.output))

    return {
        "mode" : mode,
        "inputs" : map(os.path.abspath, inputs),
        "rundir" : os.path.abspath(args.rundir),
        "reference" : os.path.abspath(args.reference),
        "noduplicates" : args.noduplicates,
        "keeptmpdir" : args.keeptmpdir,
        "barlength" : args.barcodelength,
        "bcftools" : bcftools,
        "samtools" : samtools,
        "bwa" : bwa,
        "jobs" : args.jobs,
        "output" : os.path.join(os.path.abspath(args.rundir), output)
    }
    
    
def main():
    log('Parsing and checking arguments', level='PROGRESS')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n', level='INFO')
    log('The provided input has been identified as: {}'.format(print_mode(args['mode'])), level='INFO')
    
    log('Setting up', level='PROGRESS')
    tmpdir = os.path.join(args['rundir'], '_TMP_CHISEL_PREP')
    if os.path.exists(tmpdir):
        raise ValueError("Temporary directory {} already exists, please move or rename it!".format(tmpdir))
    os.mkdir(tmpdir)
    
    if args['mode'] == 'q1':
        run = run_q1
    elif args['mode'] == 'q2':
        run = run_q2
    elif args['mode'] == 'z1':
        run = run_z1
    elif args['mode'] == 'z2':
        run = run_z2
    elif args['mode'] == 'b':
        run = run_b
    elif args['mode'] == 'B':
        run = run_B
    else:
        assert False, "Unknown mode!"
    barcoded, cells = run(args, tmpdir)
    
    log('Indexing and finalizing final barcoded BAM file', level='PROGRESS')
    indexing(args['samtools'], args['jobs'], tmpdir, barcoded, args['output'])
    
    log('Retrieving stats about the resulting barcoded BAM file', level='PROGRESS')
    cmd = '{} flagstat -@ {} {}'.format(args['samtools'], args['jobs'], args['output'])
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    log('{}'.format(stdout), level='INFO')
    
    if not args['keeptmpdir']:
        log('Cleaning remaining temporary files', level='PROGRESS')
        
    log('KTHXBYE', level='PROGRESS')
    
    
def run_q1(args, tmpdir):
    log('Running in single-end FASTQ mode', level='PROGRESS')
    
    par = {}
    par['files'] = map(lambda f : (f, ), args['inputs'])
    par['names'] = map(lambda f : os.path.basename(f[0])[:-6], par['files'])
    par['barcodes'] = mkbarcodes(par['names'], args['barlength'])
    par['tmpdir'] = tmpdir
    par['ref'] = args['reference']
    par['samtools'] = args['samtools']
    par['bwa'] = args['bwa']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Alignment, barcoding and sorting is running for every cell', level='PROGRESS')
        bams = align(**par)
    else:
        log('Alignment, barcoding, sorting, and marking duplicates is running for every cell', level='PROGRESS')
        bams = align_marking(**par)
    
    log('Merging all cells', level='PROGRESS')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir)
    
    return barcoded, list(zip(par['names'], par['barcodes']))


def run_q2(args, tmpdir):
    log('Running in paired-end FASTQ mode', level='PROGRESS')
    
    par = {}
    par['files'] = pairing_files(args['inputs'], '.fastq')
    par['names'] = map(lambda f : os.path.basename(f[0])[:-7], par['files'])
    par['barcodes'] = mkbarcodes(par['names'], args['barlength'])
    par['tmpdir'] = tmpdir
    par['ref'] = args['reference']
    par['samtools'] = args['samtools']
    par['bwa'] = args['bwa']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Alignment, barcoding and sorting is running for every cell', level='PROGRESS')
        bams = align(**par)
    else:
        log('Alignment, barcoding, sorting, and marking duplicates is running for every cell', level='PROGRESS')
        bams = align_marking(**par)
    
    log('Merging all cells', level='PROGRESS')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir)
    
    return barcoded, list(zip(par['names'], par['barcodes']))


def run_z1(args, tmpdir):
    log('Running in single-end compressed FASTQ mode', level='PROGRESS')
    
    par = {}
    par['files'] = map(lambda f : ('<(gzip -d {} )'.format(f), ), args['inputs'])
    par['names'] = map(lambda f : os.path.basename(f[0].split()[2])[:-6], par['files'])
    par['barcodes'] = mkbarcodes(par['names'], args['barlength'])
    par['tmpdir'] = tmpdir
    par['ref'] = args['reference']
    par['samtools'] = args['samtools']
    par['bwa'] = args['bwa']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Alignment, barcoding and sorting is running for every cell', level='PROGRESS')
        bams = align(**par)
    else:
        log('Alignment, barcoding, sorting, and marking duplicates is running for every cell', level='PROGRESS')
        bams = align_marking(**par)
    
    log('Merging all cells', level='PROGRESS')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir)
    
    return barcoded, list(zip(par['names'], par['barcodes']))


def run_z2(args, tmpdir):
    log('Running in paired-end compressed FASTQ mode', level='PROGRESS')
    
    par = {}
    par['files'] = map(lambda f : tuple(map(lambda p : '<(gzip -d {} )', f)), pairing_files(args['inputs'], '.fastq.gz'))
    par['names'] = map(lambda f : os.path.basename(f[0].split()[2])[:-7], par['files'])
    par['barcodes'] = mkbarcodes(par['names'], args['barlength'])
    par['tmpdir'] = tmpdir
    par['ref'] = args['reference']
    par['samtools'] = args['samtools']
    par['bwa'] = args['bwa']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Alignment, barcoding and sorting is running for every cell', level='PROGRESS')
        bams = align(**par)
    else:
        log('Alignment, barcoding, sorting, and marking duplicates is running for every cell', level='PROGRESS')
        bams = align_marked(**par)
    
    log('Merging all cells', level='PROGRESS')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir)
    
    return barcoded, list(zip(par['names'], par['barcodes']))


def run_b(args, tmpdir):
    log('Running in multiple BAM files mode', level='PROGRESS')
    
    par = {}
    par['files'] = args['inputs']
    par['names'] = map(lambda f : os.path.basename(f)[:-4], par['files'])
    par['barcodes'] = mkbarcodes(par['names'], args['barlength'])
    par['tmpdir'] = tmpdir
    par['samtools'] = args['samtools']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Barcoding and sorting is running for every cell', level='PROGRESS')
        bams = barcode(**par)
    else:
        log('Barcoding, sorting, and marking duplicates is running for every cell', level='PROGRESS')
        bams = barcode_marked(**par)
    
    log('Merging all cells', level='PROGRESS')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir)
    
    return barcoded, list(zip(par['names'], par['barcodes']))


def run_B(args, tmpdir):
    log('Running in single BAM file mode', level='PROGRESS')
    
    log('Splitting reads in BAM file by RG tag', level='PROGRESS')
    sform = os.path.join(tmpdir, '%!.bam')
    cmd = '{} split -f \'{}\' -@ {} {}'.format(args['samtools'], sform, args['jobs'], args['inputs'][0])
    stdout, stderr = sp.Popen(shlex.split(cmd), stdour=sp.PIPE, stderr=sp.PIPE).communicate()
    
    args['mode'] = 'b'
    args['input'] = map(os.path.abspath, glob.glob(os.path.join(tmpdir, '*.bam')))
    return run_b(args, tmpdir)
    
    
def mkbarcodes(names, length):
    random_sample_iter = (lambda it, k : (x for _, x in nlargest(k, ((np.random.random(), x) for x in it))))
    combr = combinations_with_replacement
    barcodes = random_sample_iter(combr(['A', 'T', 'C', 'G'], length), len(names))
    return map(lambda b : ''.join(b), barcodes)
    
    
def align(files, names, barcodes, tmpdir, ref, bwa, samtools, J):
    jobs = zip(files, names, barcodes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, bwa, ref, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_align, initargs=initargs)
    progress = (lambda e : bar.progress(advance=True, msg="Processed cell {}".format(e)))
    bams = [b for e, b in pool.imap_unordered(aligning, jobs) if progress(e)]
    pool.close()
    pool.join()
    return bams
    
    
def init_align(_tmpdir, _bwa, _ref, _samtools):
    global cmd_bwa, cmd_arg, cmd_sor, tmpdir
    cmd_bwa = '{} mem -M {} {}'.format(_bwa, _ref, '{}')
    cmd_arg = '{} addreplacerg - -r \'ID:{}\' -r \'SM:{}\''.format(_samtools, '{}', '{}')
    cmd_sor = '{} sort - -b -o {} -T {}'.format(_samtools, '{}', '{}')
    tmpdir = _tmpdir
    
    
def aligning(job):
    fil, name, barcode = job
    bam = os.path.join(tmpdir, '{}.bam'.format(name))
    curr_tmp = os.path.join(tmpdir, '_SORT_{}'.format(name))
    os.mkdir(curr_tmp)
    curr_cmd_bwa = cmd_bwa.format(' '.join(fil))
    curr_cmd_arg = cmd_arg.format('CB:Z:{}'.format(barcode), '{}_CB:Z:{}'.format(name, barcode))
    curr_cmd_sam = cmd_sor.format(bam, curr_tmp)
    pbwa = sp.Popen(shlex.split(curr_cmd_bwa), stdout=sp.PIPE, stderr=sp.PIPE)
    parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pbwa.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(shlex.split(curr_cmd_sam), stdin=parg.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return bam


def align_marked(files, names, barcodes, tmpdir, ref, bwa, samtools, J):
    jobs = zip(files, names, barcodes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, bwa, ref, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_align_marked, initargs=initargs)
    progress = (lambda e : bar.progress(advance=True, msg="Processed cell {}".format(e)))
    bams = [b for e, b in pool.imap_unordered(aligning_marked, jobs) if progress(e)]
    pool.close()
    pool.join()
    return bams
    
    
def init_align_marked(_tmpdir, _bwa, _ref, _samtools):
    global cmd_bwa, cmd_nam, cmd_fix, cmd_arg, cmd_sor, cmd_mar, tmpdir
    cmd_bwa = '{} mem -M {} {}'.format(_bwa, _ref, '{}')
    cmd_nam = '{} sort - -n -T {}'.format(_samtools, '{}')
    cmd_fix = '{} fixmate -m - -'.format(_samtools)
    cmd_arg = '{} addreplacerg - -r \'ID:{}\' -r \'SM:{}\''.format(_samtools, '{}', '{}')
    cmd_sor = '{} sort - -b -T {}'.format(_samtools, '{}')
    cmd_mar = '{} markdup -T {} - {}'.format(_samtools, '{}', '{}')
    tmpdir = _tmpdir
    
    
def aligning_marked(job):
    fil, name, barcode = job
    bam = os.path.join(tmpdir, '{}.bam'.format(name))
    nam_tmp = os.path.join(tmpdir, '_NAME_{}'.format(name))
    os.mkdir(nam_tmp)
    sor_tmp = os.path.join(tmpdir, '_SORT_{}'.format(name))
    os.mkdir(sor_tmp)
    mar_tmp = os.path.join(tmpdir, '_MARK_{}'.format(name))
    os.mkdir(mar_tmp)
    curr_cmd_bwa = cmd_bwa.format(' '.join(fil))
    curr_cmd_nam = cmd_nam.format(nam_tmp)
    curr_cmd_fix = cmd_fix
    curr_cmd_arg = cmd_arg.format('CB:Z:{}'.format(barcode), '{}_CB:Z:{}'.format(name, barcode))
    curr_cmd_sor = cmd_sor.format(sor_tmp)
    curr_cmd_mar = cmd_mar.format(mar_tmp, bam)
    pbwa = sp.Popen(shlex.split(curr_cmd_bwa), stdout=sp.PIPE, stderr=sp.PIPE)
    pnam = sp.Popen(shlex.split(curr_cmd_nam), stdin=pbwa.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    pfix = sp.Popen(shlex.split(curr_cmd_fix), srdin=pnam.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pfix.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=sp.PIPE)    
    stdout, stderr = sp.Popen(shlex.split(curr_cmd_mar), stdin=psor.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return bam


def barcode(files, names, barcodes, tmpdir, samtools, J):
    jobs = zip(files, names, barcodes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_barcoding, initargs=initargs)
    progress = (lambda e : bar.progress(advance=True, msg="Processed cell {}".format(e)))
    bams = [b for e, b in pool.imap_unordered(barcoding, jobs) if progress(e)]
    pool.close()
    pool.join()
    return bams
    
    
def init_barcoding(_tmpdir, _samtools):
    global cmd_arg, tmpdir
    cmd_arg = '{} addreplacerg {} -r \'ID:{}\' -r \'SM:{}\' -o {}'.format(_samtools, '{}', '{}', '{}', '{}')
    tmpdir = _tmpdir
    
    
def barcoding(job):
    fil, name, barcode = job
    bam = os.path.join(tmpdir, '{}.bam'.format(name))
    cmd = cmd_arg.format(fil, 'CB:Z:{}'.format(barcode), '{}_CB:Z:{}'.format(name, barcode), bam)
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return bam


def barcode_marked(files, names, barcodes, tmpdir, samtools, J):
    jobs = zip(files, names, barcodes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_barcoding_marked, initargs=initargs)
    progress = (lambda e : bar.progress(advance=True, msg="Processed cell {}".format(e)))
    bams = [b for e, b in pool.imap_unordered(barcoding_marked, jobs) if progress(e)]
    pool.close()
    pool.join()
    return bams
    
    
def init_barcoding_marked(_tmpdir, _samtools):
    global cmd_nam, cmd_fix, cmd_arg, cmd_sor, cmd_mar, tmpdir
    cmd_nam = '{} sort {} -n -T {}'.format(_samtools, '{}', '{}')
    cmd_fix = '{} fixmate -m - -'.format(_samtools)
    cmd_arg = '{} addreplacerg - -r \'ID:{}\' -r \'SM:{}\''.format(_samtools, '{}', '{}')
    cmd_sor = '{} sort - -b -T {}'.format(_samtools, '{}')
    cmd_mar = '{} markdup -T {} - {}'.format(_samtools, '{}', '{}')
    tmpdir = _tmpdir
    
    
def barcoding_marked(job):
    fil, name, barcode = job
    bam = os.path.join(tmpdir, '{}.bam'.format(name))
    nam_tmp = os.path.join(tmpdir, '_NAME_{}'.format(name))
    os.mkdir(nam_tmp)
    sor_tmp = os.path.join(tmpdir, '_SORT_{}'.format(name))
    os.mkdir(sor_tmp)
    mar_tmp = os.path.join(tmpdir, '_MARK_{}'.format(name))
    os.mkdir(mar_tmp)
    curr_cmd_nam = cmd_nam.format(fil, nam_tmp)
    curr_cmd_fix = cmd_fix
    curr_cmd_arg = cmd_arg.format('CB:Z:{}'.format(barcode), '{}_CB:Z:{}'.format(name, barcode))
    curr_cmd_sor = cmd_sor.format(sor_tmp)
    curr_cmd_mar = cmd_mar.format(mar_tmp, bam)
    pnam = sp.Popen(shlex.split(curr_cmd_nam), stdout=sp.PIPE, stderr=sp.PIPE)
    pfix = sp.Popen(shlex.split(curr_cmd_fix), srdin=pnam.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pfix.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=sp.PIPE)    
    stdout, stderr = sp.Popen(shlex.split(curr_cmd_mar), stdin=psor.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return bam


def merging(samtools, bams, jobs, tmpdir):
    barcoded = os.path.join(tmpdir, 'barcodedcells.bam')
    cellbams = os.path.join(tmpdir, 'cellbams.tsv')
    with open(cellbams, 'w') as o:
        o.write('\n'.join(bams))
    cmd = '{} merge {} -b {} -@ {}'.format(samtools, barcoded, cellbams, jobs)
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return barcoded


def indexing(samtools, jobs, tmpdir, barcoded, output):
    shutil.move(barcoded, output)
    cmd = '{} index {} -@ {}'.format(samtools, output, jobs)
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return

    
def print_mode(mode):
    if mode == 'q1':
        return 'Multiple single-ended FASTQ files!'
    elif mode == 'q2':
        return 'Multiple paired-ended FASTQ files!'
    elif mode == 'z1':
        return 'Multiple single-ended compressed FASTQ files!'
    elif mode == 'z2':
        return 'Multiple paired-ended compressed FASTQ files!'
    elif mode == 'b':
        return 'Multiple BAM files'
    elif mode == 'B':
        return 'Single BAM file'
    else:
        assert False, "Unknown mode!"
        
        
def pairing_files(files, ext):
    bases = map(lambda f : f[:-(1 + len(ext))], files)
    first = map(lambda b : [f for f in files if '{}1{}'.format(b, ext) == f], bases)
    assert all(len(f) == 1 for f in first), 'First end has not been found for {} files within:\n{}'.format(ext, files)
    second = map(lambda b : [f for f in files if '{}2{}'.format(b, ext) == f], bases)
    assert all(len(f) == 1 for f in second), 'Second end has not been found for {} files within:\n{}'.format(ext, files)
    return list(zip(first, second))
    
    
if __name__ == '__main__':
    main()
