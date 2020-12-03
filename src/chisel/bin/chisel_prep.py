#!/usr/bin/env python2.7

import os
import argparse
import shlex
import shutil
import glob
import re
import subprocess as sp
import multiprocessing as mp

from itertools import product
from collections import defaultdict
from heapq import nlargest

import numpy as np

import chisel

src = os.path.dirname(chisel.__file__)
from ..Utils import *


def parse_args():
    description = "CHISEL command to create a barcoded BAM file from single-cell FASTQs, single-cell BAMs, or a `RG:Z:`-barcoded BAM files without `CB:Z:` field."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='+', type=str, help="Input FASTQs, BAMs, or TSV file with different behaviors: ................................. (1) FASTQs -- specified in a directory DIR as `DIR/*.fastq` or `DIR/*.fastq.gz` -- will be barcoded and aligned into a barcoded BAM file where sample name, lane and read mate are extracted from filenames; .............. (2) BAMs -- specified in a directory DIR as `DIR/*.bam` -- will be barcoded and aligned into a barcoded BAM file where the cell name coincide with the name of each BAM file; ................................. (3) a single BAM file with barcodes in the field `RG:Z:` will be converted into a barcoded BAM file with the additional `CB:Z:` field. ......(4) a TSV file specifying FASTQ files (4 columns specifying FILE, SAMPLE/CELL-NAME, LANE, and READ-MATE) or BAM files (1 or 2 columns with FILE and optional SAMPLE/CELL-NAME) where header, if present, must start with symbol `#`. .....Note that cells with the same name (and different LANE for FASTQs) will be assigned the same unique cell barcode (i.e. same cell); also FASTQs with two different reads will be considered paired-end.")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-o","--output", required=False, default='barcodedcells.bam', type=str, help="Output name in running directory (default: barcodedcells.bam)")    
    parser.add_argument("--noduplicates", required=False, default=False, action='store_true', help="Do not perform marking duplicates and recalibration with Picard tools (default: False)")
    parser.add_argument("--keeptmpdir", required=False, default=False, action='store_true', help="Do not erase temporary directory (default: False)")
    parser.add_argument("--barcodelength", required=False, type=int, default=12, help="Length of barcodes (default: 12)")
    parser.add_argument("--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable (default: in $PATH)")
    parser.add_argument("--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable (default: in $PATH)")
    parser.add_argument("--bwa", required=False, default=None, type=str, help="Path to the directory to \"bwa\" executable (default: in $PATH)")
    parser.add_argument("--rexpname", required=False, default='(.*)_S.*_L.*_R[1|2]_001.fastq.*', type=str, help="Regulare expression to extract cell/sample name from input FASTQ filenames (default: `(.*)_S.*_L.*_R[1|2]_001.fastq.*`)")
    parser.add_argument("--rexplane", required=False, default='.*_S.*_(L.*)_R[1|2]_001.fastq.*', type=str, help="Regulare expression to extract cell/sample name from input FASTQ filenames (default: `.*_S.*_(L.*)_R[1|2]_001.fastq.*`)")
    parser.add_argument("--rexpread", required=False, default='.*_S.*_L.*_(R[1|2])_001.fastq.*', type=str, help="Regulare expression to extract cell/sample name from input FASTQ filenames (default: `.*_S.*_L.*_(R[1|2])_001.fastq.*`)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("--seed", required=False, type=int, default=None, help="Random seed for replication (default: None)")
    args = parser.parse_args()
    
    if args.seed is not None:
        np.random.seed(args.seed)
        
    inputs = map(os.path.abspath, filter(lambda f : len(f) > 1, args.INPUT))
    for inp in inputs:
        if not os.path.isfile(inp):
            raise ValueError("This input file does not exist: {}".format(inp))
    
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
    
    output = os.path.basename(args.output if args.output[-4:] == '.bam' else '{}.bam'.format(args.output))
    
    return {
        "inputs" : inputs,
        "rundir" : os.path.abspath(args.rundir),
        "reference" : os.path.abspath(args.reference),
        "noduplicates" : args.noduplicates,
        "keeptmpdir" : args.keeptmpdir,
        "barlength" : args.barcodelength,
        "bcftools" : bcftools,
        "samtools" : samtools,
        "bwa" : bwa,
        "rexpname" : args.rexpname,
        "rexplane" : args.rexplane,
        "rexpread" : args.rexpread,
        "jobs" : args.jobs,
        "output" : os.path.join(os.path.abspath(args.rundir), output)
    }
    
    
def main():
    log('Parsing and checking arguments', level='STEP')
    args = parse_args()
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args if a != 'inputs']) + '\n', level='INFO')
    
    log('Setting up', level='STEP')
    tmpdir = os.path.join(args['rundir'], '_TMP_CHISEL_PREP')
    if os.path.exists(tmpdir):
        raise ValueError("Temporary directory {} already exists, please move or rename it!".format(tmpdir))
    os.mkdir(tmpdir)
    errdir = os.path.join(args['rundir'], '_ERR_CHISEL_PREP')
    if os.path.exists(errdir):
        raise ValueError("Temporary error directory {} already exists, please move or rename it!".format(errdir))
    os.mkdir(errdir)
    
    if len(args['inputs']) == 1 and args['inputs'][0][:-4] == '.tsv':
        args['inputs'], info = read_table(args['inputs'][0])
    else:
        info = None

    if all(f[-6:] == '.fastq' or f[-9:] == '.fastq.gz' for f in args['inputs']):
        if info is None:
            files, fastqinfo, ispaired = match_fastq(args['inputs'], args)
        else:
            files, fastqinfo, ispaired = make_fastqinfo(args['inputs'], info, args)
        if ispaired:
            log('Running in paired-end FASTQ mode', level='STEP')
        else:
            log('Running in single-end FASTQ mode', level='STEP')
        barcoded, cells = run_q(args, tmpdir, errdir, files, fastqinfo)
        header = '#FILES\tCELL\tLANE\tREADS\tBARCODE'
        loginfo = map(lambda c : (','.join(c[0]), fastqinfo[c[0][0]]['sample'], fastqinfo[c[0][0]]['lane'], ','.join([fastqinfo[f]['read'] for f in c[0]]), c[2]), cells)
            
    elif all(f[-4:] == '.bam' for f in args['inputs']):
        if info is None:
            info = make_baminfo(args['inputs'])
        barcoded, cells = run_B(args, tmpdir, errdir) if len(args['inputs']) == 1 else run_b(args, tmpdir, errdir, binfo=info)
        header = '#FILE\tCELL\tLANE\tBARCODE'
        loginfo = map(lambda c : (c[0], info[c[0]]['sample'], info[c[0]]['lane'], c[2]), cells)  
        
    else:
        raise ValueError("Input files are of wrong format or mixed formats")
    
    log('Indexing and finalizing final barcoded BAM file', level='STEP')
    indexing(args['samtools'], args['jobs'], tmpdir, errdir, barcoded, args['output'])
    
    log('Retrieving stats about the resulting barcoded BAM file', level='STEP')
    cmd = '{} flagstat -@ {} {}'.format(args['samtools'], args['jobs'], args['output'])
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    log('{}'.format(stdout), level='INFO')
    
    floginfo = os.path.join(args['rundir'], '{}.info.tsv'.format(args['output'][-4:]))
    log('Writing final summary of barcoded cells', level='STEP')
    log('Number of barcoded cells: {}'.format(len(set(c[-1] for c in cells))), level='INFO')
    with open(floginfo, 'w') as o:
        o.write('{}\n'.format(header))
        for l in loginfo:
            o.write('{}\n'.format('\t'.join(l)))
    log('Final summary is written in {}'.format(floginfo), level='INFO')

    if not args['keeptmpdir']:
        log('Cleaning remaining temporary files', level='STEP')
        shutil.rmtree(tmpdir)
        shutil.rmtree(errdir)
        
    log('KTHXBYE', level='STEP')
    
    
def read_table(file):
    with open(file, 'r') as i:
        read = [l.strip().split() for l in i if l[0] != '#' and len(l) > 1]
    get_ext = (lambda s : 'b' if s[-4:] == '.bam' else ('q' if s[-6:] == '.fastq' or s[-9:] == '.fastq.gz' else None))
    exts = set(get_ext(r[0]) for r in read)
    if None in exts:
        raise ValueError('Unknown format, different than .bam, .fastq, or .fastq.gz has been provided!')
    elif {'b'} == exts:
        inputs = map(lambda r : os.path.abspath(r[0]), read)
        for inp in inputs:
            if not os.path.isfile(inp):
                raise ValueError("This input file does not exist: {}".format(inp))
        binfo = {os.path.abspath(r[0]) : r[1] if len(r) > 1 else os.path.abspath(r[0]) for r in read}
        lanes = defaultdict(lambda : [])
        map(lambda f : lanes[binfo[f]].append((f, len(lanes[binfo[f]]))), binfo)
        lanes = {s : dict(lanes[s]) for s in lanes}
        return inputs, {f : {'sample' : binfo[f], 'lane' : lanes[binfo[f]][f]} for f in binfo}
    elif {'q'} == exts:
        wrong = filter(lambda r : len(r) < 4, read)
        if len(wrong) > 0:
            raise ValueError('When FASTQ files are provided in a table, each row should have four tab-separated entries: FILE, SAMPLE/CELL-NAME, LANE, and READ!')
        inputs = map(lambda r : os.path.abspath(r[0]), read)
        for inp in inputs:
            if not os.path.isfile(inp):
                raise ValueError("This input file does not exist: {}".format(inp))
        return inputs, {os.path.abspath(r[0]) : {'sample' : r[1], 'lane' : r[2], 'read' : r[3]} for r in read}
    else:
        raise ValueError('The provided table contains files of mixed format (i.e. both BAM and FASTQ)!')


def match_fastq(inputs, args):
    mname = (lambda s : re.search(args['rexpname'], os.path.basename(s)))
    mlane = (lambda s : re.search(args['rexplane'], os.path.basename(s)))
    mread = (lambda s : re.search(args['rexpread'], os.path.basename(s)))
    match = map(lambda f : (f, mname(f), mlane(f), mread(f)), inputs)
    if all(None not in m for m in match):
        fastqinfo = {m[0] : {'sample' : m[1].group(1), 'lane' : m[2].group(1), 'read' : m[3].group(1)} for m in match}
        return make_fastqinfo(inputs, fastqinfo, args)
    else:
        log('The filesnames of the provided FASTQ files do not match the format with given expressions and will be all considered as independent single-end FASTQs', level='WARN')
        rmext = (lambda s : s[-9:] if s[-9:] == '.fastq.gz' else s[-6:])
        fastqinfo = {m[0] : {'sample' : rmext(m[0]), 'lane' : 'L001', 'read' : 'R1'} for m in match}
        files = map(lambda f : (f, ), inputs)
        return files, fastqinfo, False
    
        
def make_fastqinfo(inputs, fastqinfo, args):
    pairs = defaultdict(lambda : [])
    map(lambda f : pairs[(fastqinfo[f]['sample'], fastqinfo[f]['lane'])].append(f), fastqinfo)
    map(lambda p : pairs[p].sort(key=(lambda f : fastqinfo[f]['read'])), pairs)
    if all(len(pairs[p]) == 1 for p in pairs):
        files = map(lambda f : (f, ), inputs)
        fastqinfo = {f : fastqinfo[f] for f in fastqinfo}
        return files, fastqinfo, False
    elif all(len(set(pairs[p])) == 2 and len(pairs[p]) == 2 for p in pairs):
        files = map(lambda p : (pairs[p][0], pairs[p][1]), pairs)
        fastqinfo = {f : fastqinfo[f] for f in fastqinfo}
        assert set(f for p in files for f in p) == set(fastqinfo.keys())
        return files, fastqinfo, True
    else:
        for p in pairs:
            if set(len(pairs[p])) < 2:
                raise ValueError('Found more than TWO files with the same sample and lane but also identical reads!\n{}'.format(','.join(pairs[p])))
            elif len(pairs[p]) > 2:
                raise ValueError('Found more than TWO files with the same sample and lane, which cannot indicate paired-end reads!\n{}'.format(','.join(pairs[p])))
            else:
                assert False
                
                
def make_baminfo(inputs):
    binfo = {os.path.abspath(f) : os.path.basename(f) for f in inputs}
    lanes = defaultdict(lambda : [])
    map(lambda f : lanes[binfo[f]].append((f, len(lanes[binfo[f]]))), binfo)
    lanes = {s : dict(lanes[s]) for s in lanes}
    return {f : {'sample' : binfo[f], 'lane' : lanes[binfo[f]][f]} for f in binfo}
    

def run_q(args, tmpdir, errdir, files, fastqinfo):
    par = {}
    par['files'] = files
    par['names'] = map(lambda f : fastqinfo[f[0]]['sample'], files)
    par['lanes'] = map(lambda f : fastqinfo[f[0]]['lane'], files)
    par['barcodes'] = mkbarcodes(par['files'], args['barlength'], qinfo=fastqinfo)
    par['tmpdir'] = tmpdir
    par['errdir'] = errdir
    par['ref'] = args['reference']
    par['samtools'] = args['samtools']
    par['bwa'] = args['bwa']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Alignment, barcoding and sorting is running for every cell', level='STEP')
        bams = align(**par)
    else:
        log('Alignment, barcoding, sorting, and marking duplicates is running for every cell', level='STEP')
        bams = align_marked(**par)
    
    log('Merging all cells', level='STEP')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir, errdir)
    
    return barcoded, list(zip(par['files'], par['names'], par['barcodes']))


def run_b(args, tmpdir, errdir, binfo):
    log('Running in multiple BAM files mode', level='STEP')
    
    par = {}
    par['files'] = args['inputs']
    par['names'] = map(lambda f : binfo[f]['sample'], par['files'])
    par['lanes'] = map(lambda f : binfo[f]['lane'], par['files'])
    par['barcodes'] = mkbarcodes(par['files'], args['barlength'], binfo=binfo)
    par['tmpdir'] = tmpdir
    par['errdir'] = errdir
    par['samtools'] = args['samtools']
    par['J'] = args['jobs']
    if args['noduplicates']:
        log('Barcoding and sorting is running for every cell', level='STEP')
        bams = barcode(**par)
    else:
        log('Barcoding, sorting, and marking duplicates is running for every cell', level='STEP')
        bams = barcode_marked(**par)
    
    log('Merging all cells', level='STEP')
    barcoded = merging(args['samtools'], bams, args['jobs'], tmpdir, errdir)
    
    return barcoded, list(zip(par['files'], par['names'], par['barcodes']))


def run_B(args, tmpdir, errdir):
    log('Running in single BAM file mode', level='STEP')
    
    log('Splitting reads in BAM file by RG tag', level='STEP')
    sform = os.path.join(tmpdir, '%!.bam')
    cmd = '{} split -f \'{}\' -@ {} {}'.format(args['samtools'], sform, args['jobs'], args['inputs'][0])
    stdout, stderr = sp.Popen(shlex.split(cmd), stdour=sp.PIPE, stderr=sp.PIPE).communicate()
    
    args['input'] = map(os.path.abspath, glob.glob(os.path.join(tmpdir, '*.bam')))
    getname = (lambda f : os.path.splitext(os.path.basename(f))[0])
    binfo = {f : {'sample' : getname(f), 'lane' : 'L001'} for f in args['input']}
    return run_b(args, tmpdir, errdir, binfo)
    
    
def mkbarcodes(files, length, qinfo=None, binfo=None):
    random_sample_iter = (lambda it, k : (x for _, x in nlargest(k, ((np.random.random(), x) for x in it))))
    barcodes = random_sample_iter(product(['A', 'T', 'C', 'G'], repeat=length), len(files))
    barcodes = map(lambda b : ''.join(b), barcodes)
    assert len(barcodes) == len(files), '{} != {}'.format(len(files), len(barcodes))
    if qinfo is not None:
        lanes = defaultdict(lambda : [])
        map(lambda f : lanes[qinfo[f[0]]['sample']].append(f), files)
        dup = [s for s in lanes if len(lanes[s]) != len(set(qinfo[f[0]]['lane'] for f in lanes[s]))]
        if len(dup) > 0:
            raise ValueError('Two or more of these files have the same sample name but also the same lane number:\n{}'.format('\n'.join([f for f in lanes[dup[0]]])))
        assign = dict(zip(files, barcodes))
        barcodes = map(lambda f : assign[lanes[qinfo[f[0]]['sample']][0]], files)
    if binfo is not None:
        lanes = defaultdict(lambda : [])
        map(lambda f : lanes[binfo[f]].append(f), files)
        assign = dict(zip(files, barcodes))
        barcodes = map(lambda f : assign[lanes[binfo[f]][0]], files)        
    return barcodes
    
    
def align(files, names, barcodes, lanes, tmpdir, errdir, ref, bwa, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, bwa, ref, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_align, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c==0 else error(e, errdir))
    bams = [b for e, b, c in pool.imap_unordered(aligning, jobs) if progress(e, c)]
    pool.close()
    pool.join()
    return bams
    
    
def init_align(_tmpdir, _errdir, _bwa, _ref, _samtools):
    global cmd_bwa, cmd_arg, cmd_sor, tmpdir, errdir
    cmd_bwa = '{} mem -M {} {}'.format(_bwa, _ref, '{}')
    cmd_arg = '{} addreplacerg - -r \'ID:{}\' -r \'SM:{}\' -r \'FO:{}\' -r \'PG:CHISEL_PREP\' -Osam'.format(_samtools, '{}', '{}', '{}')
    cmd_sor = '{} sort - -Obam -o {} -T {}'.format(_samtools, '{}', '{}')
    tmpdir = _tmpdir
    errdir = _errdir
    
    
def aligning(job):
    fil, name, barcode, lane = job
    bam = os.path.join(tmpdir, '{}_{}.bam'.format(name, lane))
    curr_tmp = os.path.join(tmpdir, '_SORT_{}_{}'.format(name, lane))
    os.mkdir(curr_tmp)
    curr_cmd_bwa = cmd_bwa.format(' '.join(fil))
    curr_cmd_arg = cmd_arg.format('CB:Z:{}-{}'.format(barcode, lane), '{}'.format(name), '{}'.format(lane))
    curr_cmd_sam = cmd_sor.format(bam, curr_tmp)
    pbwa = sp.Popen(shlex.split(curr_cmd_bwa), stdout=sp.PIPE, stderr=sp.PIPE)
    parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pbwa.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(shlex.split(curr_cmd_sam), stdin=parg.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    with open(os.path.join(errdir, '{}_{}.log'.format(name, lane))) as o:
        o.write('#### BWA ####\n{}\n'.format(pbwa.stderr))
        o.write('#### SAMTOOLS ADDREPLACERG ####\n{}\n'.format(parg.stderr))
        o.write('#### SAMTOOLS SORT ####\n{}\n'.format(stderr))
    return '{}_{}'.format(name, lane), bam


def align_marked(files, names, barcodes, lanes, tmpdir, errdir, ref, bwa, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, bwa, ref, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_align_marked, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c==0 else error(e, errdir))
    bams = [b for e, b, c in pool.imap_unordered(aligning, jobs) if progress(e, c)]
    pool.close()
    pool.join()
    return bams
    
    
def init_align_marked(_tmpdir, _errdir, _bwa, _ref, _samtools):
    global cmd_bwa, cmd_nam, cmd_fix, cmd_arg, cmd_sor, cmd_mar, tmpdir, errdir
    cmd_bwa = '{} mem -M {} {}'.format(_bwa, _ref, '{}')
    cmd_nam = '{} sort - -n -T {} -Osam'.format(_samtools, '{}')
    cmd_fix = '{} fixmate -m - - -Osam'.format(_samtools)
    cmd_arg = '{} addreplacerg - -r \'ID:{}\' -r \'SM:{}\' -r \'FO:{}\' -r \'PG:CHISEL_PREP\' -Osam'.format(_samtools, '{}', '{}', '{}')
    cmd_sor = '{} sort - -T {} -Osam'.format(_samtools, '{}')
    cmd_mar = '{} markdup -T {} - {} -Obam'.format(_samtools, '{}', '{}')
    tmpdir = _tmpdir
    errdir = _errdir
    
    
def aligning_marked(job):
    fil, name, barcode, lane = job
    bam = os.path.join(tmpdir, '{}_{}.bam'.format(name, lane))
    nam_tmp = os.path.join(tmpdir, '_NAME_{}_{}'.format(name, lane))
    os.mkdir(nam_tmp)
    sor_tmp = os.path.join(tmpdir, '_SORT_{}_{}'.format(name, lane))
    os.mkdir(sor_tmp)
    mar_tmp = os.path.join(tmpdir, '_MARK_{}_{}'.format(name, lane))
    os.mkdir(mar_tmp)
    curr_cmd_bwa = cmd_bwa.format(' '.join(fil))
    curr_cmd_nam = cmd_nam.format(nam_tmp)
    curr_cmd_fix = cmd_fix
    curr_cmd_arg = cmd_arg.format('CB:Z:{}-{}'.format(barcode, lane), '{}'.format(name), '{}'.format(lane))
    curr_cmd_sor = cmd_sor.format(sor_tmp)
    curr_cmd_mar = cmd_mar.format(mar_tmp, bam)
    pbwa = sp.Popen(shlex.split(curr_cmd_bwa), stdout=sp.PIPE, stderr=sp.PIPE)
    pnam = sp.Popen(shlex.split(curr_cmd_nam), stdin=pbwa.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    pfix = sp.Popen(shlex.split(curr_cmd_fix), stdin=pnam.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pfix.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=sp.PIPE)    
    stdout, stderr = sp.Popen(shlex.split(curr_cmd_mar), stdin=psor.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    with open(os.path.join(errdir, '{}_{}.log'.format(name, lane))) as o:
        o.write('#### BWA ####\n{}\n'.format(pbwa.stderr))
        o.write('#### SAMTOOLS SORT BY NAME ####\n{}\n'.format(pnam.stderr))
        o.write('#### SAMTOOLS FIXMATE ####\n{}\n'.format(pfix.stderr))
        o.write('#### SAMTOOLS ADDREPLACERG ####\n{}\n'.format(parg.stderr))
        o.write('#### SAMTOOLS SORT ####\n{}\n'.format(psor.stderr))
        o.write('#### SAMTOOLS MARKDUP ####\n{}\n{}\n'.format(stderr, stdout))
    return '{}_{}'.format(name, lane), bam


def barcode(files, names, barcodes, lanes, tmpdir, errdir, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_barcoding, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c==0 else error(e, errdir))
    bams = [b for e, b, c in pool.imap_unordered(aligning, jobs) if progress(e, c)]
    pool.close()
    pool.join()
    return bams
    
    
def init_barcoding(_tmpdir, _errdir, _samtools):
    global cmd_arg, tmpdir, errdir
    cmd_arg = '{} addreplacerg {} -r \'ID:{}\' -r \'SM:{}\' -r \'FO:{}\' -r \'PG:CHISEL_PREP\' -o {}'.format(_samtools, '{}', '{}', '{}', '{}', '{}')
    tmpdir = _tmpdir
    errdir = _errdir
    
    
def barcoding(job):
    fil, name, barcode, lane = job
    bam = os.path.join(tmpdir, '{}_{}.bam'.format(name, lane))
    cmd = cmd_arg.format(fil, 'CB:Z:{}-{}'.format(barcode, lane), '{}'.format(name), '{}'.format(lane), bam)
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    with open(os.path.join(errdir, '{}_{}.log'.format(name, lane))) as o:
        o.write('#### SAMTOOLS ADDREPLACERG ####\n{}\n{}\n'.format(stderr, stdout))
    return name, bam


def barcode_marked(files, names, barcodes, lanes, tmpdir, errdir, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_barcoding_marked, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c==0 else error(e, errdir))
    bams = [b for e, b, c in pool.imap_unordered(aligning, jobs) if progress(e, c)]
    pool.close()
    pool.join()
    return bams
    
    
def init_barcoding_marked(_tmpdir, _errdir, _samtools):
    global cmd_nam, cmd_fix, cmd_arg, cmd_sor, cmd_mar, tmpdir, errdir
    cmd_nam = '{} sort {} -n -T {} -Osam'.format(_samtools, '{}', '{}')
    cmd_fix = '{} fixmate -m - - -Osam'.format(_samtools)
    cmd_arg = '{} addreplacerg - -r \'ID:{}\' -r \'SM:{}\' -r \'FO:{}\' -r \'PG:CHISEL_PREP\' -Osam'.format(_samtools, '{}', '{}', '{}')
    cmd_sor = '{} sort - -T {} -Osam'.format(_samtools, '{}')
    cmd_mar = '{} markdup -T {} - {} -Obam'.format(_samtools, '{}', '{}')
    tmpdir = _tmpdir
    errdir = _errdir
    
    
def barcoding_marked(job):
    fil, name, barcode, lane = job
    bam = os.path.join(tmpdir, '{}_{}.bam'.format(name, lane))
    nam_tmp = os.path.join(tmpdir, '_NAME_{}_{}'.format(name, lane))
    os.mkdir(nam_tmp)
    sor_tmp = os.path.join(tmpdir, '_SORT_{}_{}'.format(name, lane))
    os.mkdir(sor_tmp)
    mar_tmp = os.path.join(tmpdir, '_MARK_{}_{}'.format(name, lane))
    os.mkdir(mar_tmp)
    curr_cmd_nam = cmd_nam.format(fil, nam_tmp)
    curr_cmd_fix = cmd_fix
    curr_cmd_arg = cmd_arg.format('CB:Z:{}-{}'.format(barcode, lane), '{}'.format(name), '{}'.format(lane))
    curr_cmd_sor = cmd_sor.format(sor_tmp)
    curr_cmd_mar = cmd_mar.format(mar_tmp, bam)
    pnam = sp.Popen(shlex.split(curr_cmd_nam), stdout=sp.PIPE, stderr=sp.PIPE)
    pfix = sp.Popen(shlex.split(curr_cmd_fix), stdin=pnam.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pfix.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=sp.PIPE)    
    stdout, stderr = sp.Popen(shlex.split(curr_cmd_mar), stdin=psor.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    with open(os.path.join(errdir, '{}_{}.log'.format(name, lane))) as o:
        o.write('#### SAMTOOLS SORT BY NAME ####\n{}\n'.format(pnam.stderr))
        o.write('#### SAMTOOLS FIXMATE ####\n{}\n'.format(pfix.stderr))
        o.write('#### SAMTOOLS ADDREPLACERG ####\n{}\n'.format(parg.stderr))
        o.write('#### SAMTOOLS SORT ####\n{}\n'.format(psor.stderr))
        o.write('#### SAMTOOLS MARKDUP ####\n{}\n{}\n'.format(stderr, stdout))
    return name, bam


def merging(samtools, bams, jobs, tmpdir, errdir):
    barcoded = os.path.join(tmpdir, 'barcodedcells.bam')
    cellbams = os.path.join(tmpdir, 'cellbams.tsv')
    with open(cellbams, 'w') as o:
        o.write('\n'.join(bams))
    cmd = '{} merge {} -b {} -@ {}'.format(samtools, barcoded, cellbams, jobs)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise ValueError('Merging failed with messages:\n{}\n{}\n'.format(stdout, stderr))
    return barcoded


def indexing(samtools, jobs, tmpdir, errdir, barcoded, output):
    shutil.move(barcoded, output)
    cmd = '{} index {} -@ {}'.format(samtools, output, jobs)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise ValueError('Indexing failed with messages:\n{}\n{}\n'.format(stdout, stderr))
    return


def error_code(e, errdir):
    raise ValueError('Commands failed for {}, check the errors in the corresponding log file in:\n{}'.format(e, errdir))

    
if __name__ == '__main__':
    main()
