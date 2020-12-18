#!/usr/bin/env python2.7

import sys, os
import argparse
import shlex
import shutil
import glob
import re
import traceback
import subprocess as sp
import multiprocessing as mp

from itertools import product
from collections import defaultdict
from heapq import nlargest

import numpy as np

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'src')
if not os.path.isdir(src):
    raise ValueError(error("src directory not found in parent directory of bin i.e. {}, is anything been moved?".format(src)))
sys.path.append(src)
from Utils import *


def parse_args():
    description = '''CHISEL command to create a barcoded BAM file from single-cell FASTQs (or gz-compressed FASTQs), single-cell BAMs, or a
        `RG:Z:`-barcoded BAM files without `CB:Z:` tags. When single-cell FASTQs or BAMs are provided
        a CELL name is assigned to each file (through either filename or table) and the same cell barcode will be assigned to all corresponding reads, but
        a different RG tag as they are considered as different repetitions of sequencing of the same cell. Specifically, when a table of inputs is not provied,
        for FASTQs each CELL name is extracted from the filename through the provided regular expression (default matches Illumina standard format), for BAMs
        basename is used as CELL name. When single-cell FASTQs are provided a READ value is also assigned to each file (through either filename or table) and 
        files with the same filename when removing READ values are considered as pairs of sequencing read mates.
        Input files, CELL names, and possible READ values can be provided through a table of inputs.'''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", nargs='+', type=str, help='''Input FASTQs, BAMs, or TSV file with different behaviors: ......................................... 
    (1) FASTQs -- specified in a directory DIR as `DIR/*.fastq` or `DIR/*.fastq.gz` -- will be barcoded and aligned with (optionally) marked duplicates into a barcoded BAM file;
    .................................
    (2) BAMs -- specified in a directory DIR as `DIR/*.bam` -- will be barcoded and aligned with (optionally) marked duplicates into a barcoded BAM file; 
    ..............................................
    (3) a single BAM file with unique cells names in the field `RG:Z:` will be converted into a barcoded BAM file with the additional `CB:Z:` tag; 
    ..............
    (4) a tab-separated table of inputs (TSV with optional header starting with `#`) with two columns: the first column is an input file (FASTQ or BAM) and the
    second column is the corresponding cell name. When FASTQs are provided, a third column can be optionally specified to indicate the read name in paired-end
    sequencing, e.g., indicating either R1 or R2 for the first or second mate of paired-end reads, respectively. If a third column is not present, FASTQs are
    assumed to be from single-end sequencing.''')
    parser.add_argument("-r","--reference", type=str, required=False, help="Reference genome, which is mandatory in FASTQ mode (default: None)")
    parser.add_argument("-x","--rundir", required=False, default='./', type=str, help="Running directory (default: current directory)")
    parser.add_argument("-o","--output", required=False, default='barcodedcells.bam', type=str, help="Output name in running directory (default: barcodedcells.bam)")    
    parser.add_argument("--rexpname", required=False, default='(.*)_S.*_L00.*_R[1|2]_001.fastq.*', type=str, help="Regulare expression to extract cell name from input FASTQ filenames (default: `(.*)_S.*_L.*_R[1|2]_001.fastq.*`)")
    parser.add_argument("--rexpread", required=False, default='.*_S.*_L00.*_(R[1|2])_001.fastq.*', type=str, help="Regulare expression to extract cell name from input FASTQ filenames (default: `.*_S.*_L.*_(R[1|2])_001.fastq.*`)")
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
        
    inputs = map(os.path.abspath, filter(lambda f : len(f) > 1, args.INPUT))
    for inp in inputs:
        if not (os.path.exists(inp) and not os.path.isdir(inp)):
            raise ValueError(error("This input file does not exist: {}".format(inp)))
    
    if not os.path.isdir(args.rundir):
        raise ValueError(error("Running directory does not exists: {}".format(args.rundir)))
    if args.reference is not None:
        if not os.path.isfile(args.reference):
            raise ValueError(error("Reference genome file does not exist: {}".format(args.reference)))
        refidx = ['{}.{}'.format(args.reference, ix) for ix in ['amb', 'ann', 'bwt', 'pac', 'sa']]
        if not all(os.path.isfile(f) for f in refidx):
            raise ValueError(error("Some of the BWA index files are missing, please make sure these are available and generated through the command ``bwa index FASTA-REFERENCE''. Expected files are: {}".format('\n'.join(refidx))))

    if not args.jobs:
        args.jobs = mp.cpu_count()
    if args.jobs < 1:
        raise ValueError(error("The number of jobs must be positive!"))

    bcftools = args.bcftools
    if not bcftools:
        bcftools = "bcftools"
    if which(bcftools) is None:
        raise ValueError(error("bcftools has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} bcftools\n\nOtherwise, please provide with the flag `--bcftools` the full path to the directory containing bcftools exacutable."))

    samtools = args.samtools
    if not samtools:
        samtools = "samtools"
    if which(samtools) is None:
        raise ValueError(error("samtools has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} samtools\n\nOtherwise, please provide with the flag `--samtools` the full path to the directory containing samtools exacutable."))

    bwa = args.bwa
    if not bwa:
        bwa = "bwa"
    if which(bwa) is None:
        raise ValueError(error("bwa has not been found or is not executable!\n\nIf you are within a CHISEL conda environment ${ENV} you can install it with:\n\tconda install -c bioconda -n ${ENV} bwa\n\nOtherwise, please provide with the flag `--bwa` the full path to the directory containing bwa exacutable."))
    
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
        raise ValueError(error("Temporary directory {} already exists, please move or rename it!".format(tmpdir)))
    os.mkdir(tmpdir)
    errdir = os.path.join(args['rundir'], '_ERR_CHISEL_PREP')
    if os.path.exists(errdir):
        raise ValueError(error("Temporary error directory {} already exists, please move or rename it!".format(errdir)))
    os.mkdir(errdir)
    
    if len(args['inputs']) == 1 and args['inputs'][0][-4:] != '.bam'\
                                and args['inputs'][0][-6:] != '.fastq'\
                                and args['inputs'][0][-9:] != '.fastq.gz'\
                                and args['inputs'][0][-3:] != '.gz':
        log('Reading provided table of inputs', level='STEP')
        args['inputs'], info = read_table(args['inputs'][0])
    else:
        info = None

    if all(f[-6:] == '.fastq' or f[-9:] == '.fastq.gz' for f in args['inputs']):
        if args['reference'] is None:
            shutil.rmtree(tmpdir)
            shutil.rmtree(errdir)
            raise ValueError(error('Reference genome is required when running in FASTQ mode!'))
        if info is None:
            files, fastqinfo, ispaired = match_fastq(args['inputs'], args)
        else:
            files, fastqinfo, ispaired = make_fastqinfo(args['inputs'], info, args)
        if ispaired:
            log('Running in paired-end FASTQ mode', level='STEP')
        else:
            log('Running in single-end FASTQ mode', level='STEP')
        barcoded, cells = run_q(args, tmpdir, errdir, files, fastqinfo)
        header = '#CELL\tBARCODE\tREPETITION\tREADS\tFILES'
        cells = sorted(cells, key=(lambda c : (c[2], fastqinfo[c[0][0]]['reps'])))
        loginfo = map(lambda c : (fastqinfo[c[0][0]]['cell'], c[2], fastqinfo[c[0][0]]['reps'], ','.join([fastqinfo[f]['read'] for f in c[0]]), ','.join(c[0])), cells)
            
    elif all(f[-4:] == '.bam' for f in args['inputs']):
        if info is None:
            info = make_baminfo(args['inputs'])
        barcoded, cells, info = run_B(args, tmpdir, errdir) if len(args['inputs']) == 1 else run_b(args, tmpdir, errdir, info)
        header = '#CELL\tBARCODE\tREPETITION\tFILE'
        cells = sorted(cells, key=(lambda c : (c[2], info[c[0]]['reps'])))
        loginfo = map(lambda c : (info[c[0]]['cell'], c[2], info[c[0]]['reps'], c[0]), cells)
        
    else:
        raise ValueError(error("Input files are of wrong format or mixed formats"))
    
    log('Indexing and finalizing final barcoded BAM file', level='STEP')
    indexing(args['samtools'], args['jobs'], tmpdir, errdir, barcoded, args['output'])
    
    log('Retrieving stats about the resulting barcoded BAM file', level='STEP')
    cmd = '{} flagstat -@ {} {}'.format(args['samtools'], args['jobs'], args['output'])
    stdout, stderr = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    log('{}'.format(stdout), level='INFO')
    
    floginfo = os.path.join(args['rundir'], '{}.info.tsv'.format(args['output'][:-4]))
    log('Writing final summary of barcoded cells', level='STEP')
    log('Number of barcoded cells: {}'.format(len(set(c[-1] for c in cells))), level='INFO')
    with open(floginfo, 'w') as o:
        o.write('{}\n'.format(header))
        for l in loginfo:
            o.write('{}\n'.format('\t'.join(map(str, l))))
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
        raise ValueError(error('Unknown format, different than .bam, .fastq, or .fastq.gz has been provided!'))
    wrong = filter(lambda r : len(r) < 2, read)
    if len(wrong) > 0:
        raise ValueError(error('When a TSV is provided, every row must have at least two fields, specifying FILE and CELL-NAME, but an error was found, e.g.:\n\n{}'.format(wrong[0])))
    inputs = map(lambda r : os.path.abspath(r[0]), read)
    for inp in inputs:
        if not os.path.isfile(inp):
            raise ValueError(error("This input file does not exist: {}".format(inp)))
    if {'b'} == exts:
        binfo = {os.path.abspath(r[0]) : r[1] for r in read}
        reps = defaultdict(lambda : [])
        map(lambda f : reps[binfo[f]].append((f, len(reps[binfo[f]]))), binfo)
        reps = {s : dict(reps[s]) for s in reps}
        return inputs, {f : {'cell' : binfo[f], 'reps' : reps[binfo[f]][f]} for f in binfo}
    elif {'q'} == exts:
        qinfo = {os.path.abspath(r[0]) : {'cell' : r[1], 'read' : r[2] if len(r) > 2 else 'R1'} for r in read}
        reps = defaultdict(lambda : [])
        found = set()
        proc = (lambda Q, f, R : reps[Q['cell']].append((R, len(reps[Q['cell']]))) or found.add(R) if R not in found else None)
        map(lambda f : proc(qinfo[f], f, f.replace(qinfo[f]['read'], '')), qinfo)
        reps = {s : dict(reps[s]) for s in reps}
        return inputs, {f : {'cell' : qinfo[f]['cell'], 'reps' : reps[qinfo[f]['cell']][f.replace(qinfo[f]['read'], '')], 'read' : qinfo[f]['read']} for f in qinfo}
    else:
        raise ValueError(error('The provided table contains files of mixed format (i.e. both BAM and FASTQ)!'))


def match_fastq(inputs, args):
    mname = (lambda s : re.search(args['rexpname'], os.path.basename(s)))
    mread = (lambda s : re.search(args['rexpread'], os.path.basename(s)))
    match = map(lambda f : (f, mname(f), mread(f)), inputs)
    if all(None not in m for m in match):
        qinfo = {m[0] : {'cell' : m[1].group(1), 'read' : m[2].group(1)} for m in match}
        reps = defaultdict(lambda : [])
        found = set()
        proc = (lambda Q, f, R : reps[Q['cell']].append((R, len(reps[Q['cell']]))) or found.add(R) if R not in found else None)
        map(lambda f : proc(qinfo[f], f, f.replace(qinfo[f]['read'], '')), qinfo)
        reps = {s : dict(reps[s]) for s in reps}
        qinfo = {f : {'cell' : qinfo[f]['cell'], 'reps' : reps[qinfo[f]['cell']][f.replace(qinfo[f]['read'], '')], 'read' : qinfo[f]['read']} for f in qinfo}
        return make_fastqinfo(inputs, qinfo, args)
    else:
        log('The filesnames of the provided FASTQ files do not match the format with given expressions and will be all considered as independent single-end FASTQs', level='WARN')
        rmext = (lambda s : s[-9:] if s[-9:] == '.fastq.gz' else s[-6:])
        qinfo = {m[0] : {'cell' : rmext(m[0]), 'reps' : 0, 'read' : 'R1'} for m in match}
        files = map(lambda f : (f, ), inputs)
        return files, qinfo, False
    
        
def make_fastqinfo(inputs, fastqinfo, args):
    pairs = defaultdict(lambda : [])
    map(lambda f : pairs[(fastqinfo[f]['cell'], fastqinfo[f]['reps'])].append(f), fastqinfo)
    map(lambda p : pairs[p].sort(key=(lambda f : fastqinfo[f]['read'])), pairs)
    if all(len(pairs[p]) == 1 for p in pairs):
        files = map(lambda f : (f, ), inputs)
        return files, fastqinfo, False
    elif all(len(set(pairs[p])) == 2 and len(pairs[p]) == 2 for p in pairs):
        files = map(lambda p : (pairs[p][0], pairs[p][1]), pairs)
        assert set(f for p in files for f in p) == set(fastqinfo.keys())
        return files, fastqinfo, True
    else:
        for p in pairs:
            if set(len(pairs[p])) < 2:
                raise ValueError(error('Found more files from the same cell with same filenames and identical reads!\n{}'.format(','.join(pairs[p]))))
            elif len(pairs[p]) > 2:
                raise ValueError(error('Found more than TWO files with the same cell and filenames, which cannot indicate paired-end reads!\n{}'.format(','.join(pairs[p]))))
            else:
                assert False


def make_baminfo(inputs):
    binfo = {os.path.abspath(f) : os.path.basename(f) for f in inputs}
    lanes = defaultdict(lambda : [])
    map(lambda f : lanes[binfo[f]].append((f, len(lanes[binfo[f]]))), binfo)
    lanes = {s : dict(lanes[s]) for s in lanes}
    return {f : {'cell' : binfo[f], 'reps' : lanes[binfo[f]][f]} for f in binfo}


def run_q(args, tmpdir, errdir, files, fastqinfo):
    par = {}
    par['files'] = files
    par['names'] = map(lambda f : fastqinfo[f[0]]['cell'], files)
    par['lanes'] = map(lambda f : fastqinfo[f[0]]['reps'], files)
    par['barcodes'] = mkbarcodes(par['files'], args['barlength'], fastqinfo)
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
    par['names'] = map(lambda f : binfo[f]['cell'], par['files'])
    par['lanes'] = map(lambda f : binfo[f]['reps'], par['files'])
    par['barcodes'] = mkbarcodes(par['files'], args['barlength'], binfo)
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
    
    return barcoded, list(zip(par['files'], par['names'], par['barcodes'])), binfo


def run_B(args, tmpdir, errdir):
    log('Running in single BAM file mode', level='STEP')
    
    log('Splitting reads in BAM file by RG tag', level='STEP')
    sform = os.path.join(tmpdir, '%!.bam')
    cmd = '{} split -f \'{}\' -@ {} {}'.format(args['samtools'], sform, args['jobs'], args['inputs'][0])
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise ValueError(error('Merging failed with messages:\n{}\n{}\n'.format(stdout, stderr)))
    
    args['inputs'] = map(os.path.abspath, glob.glob(os.path.join(tmpdir, '*.bam')))
    getname = (lambda f : os.path.splitext(os.path.basename(f))[0])
    binfo = {f : {'cell' : getname(f), 'reps' : 0} for f in args['inputs']}
    return run_b(args, tmpdir, errdir, binfo)
    
    
def mkbarcodes(files, length, info):
    random_sample_iter = (lambda it, k : (x for _, x in nlargest(k, ((np.random.random(), x) for x in it))))
    barcodes = random_sample_iter(product(['A', 'T', 'C', 'G'], repeat=length), len(files))
    barcodes = map(lambda b : ''.join(b), barcodes)
    assert len(barcodes) == len(files), '{} != {}'.format(len(files), len(barcodes))
    getfile = (lambda f : f[0] if type(f) == tuple else f)
    reps = defaultdict(lambda : [])
    map(lambda f : reps[info[getfile(f)]['cell']].append(f), files)
    dup = [s for s in reps if len(reps[s]) != len(set(info[getfile(f)]['reps'] for f in reps[s]))]
    if len(dup) > 0:
        raise ValueError(error('Two or more of these files have the same cell name but also the same rep number:\n{}'.format('\n'.join([f for f in reps[dup[0]]]))))
    assign = dict(zip(files, barcodes))
    return map(lambda f : assign[reps[info[getfile(f)]['cell']][0]], files)
    
    
def align(files, names, barcodes, lanes, tmpdir, errdir, ref, bwa, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, bwa, ref, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_align, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c is None else error_code(e, c, errdir))
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
    curr_cmd_sor = cmd_sor.format(bam, curr_tmp)
    blog = os.path.join(errdir, '{}_{}_BWA.log'.format(name, lane))
    rlog = os.path.join(errdir, '{}_{}_SAMTOOLS_ADDREPLACERG.log'.format(name, lane))
    olog = os.path.join(errdir, '{}_{}_SAMTOOLS_SORT.log'.format(name, lane))
    with open(blog, 'w') as ebwa, open(rlog, 'w') as earg, open(olog, 'w') as osor:
        pbwa = sp.Popen(shlex.split(curr_cmd_bwa), stdout=sp.PIPE, stderr=ebwa)
        parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pbwa.stdout, stdout=sp.PIPE, stderr=earg)
        psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=osor)
        rcodes = map(lambda p : p.wait(), [pbwa, parg, psor])
    return '{}_{}'.format(name, lane), bam, check_rcodes(rcodes, [blog, rlog, olog], [curr_tmp])


def align_marked(files, names, barcodes, lanes, tmpdir, errdir, ref, bwa, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, bwa, ref, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_align_marked, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c is None else error_code(e, c, errdir))
    bams = [b for e, b, c in pool.imap_unordered(aligning_marked, jobs) if progress(e, c)]
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
    cmd_mar = '{} markdup -T {} -Obam - {}'.format(_samtools, '{}', '{}')
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
    blog = os.path.join(errdir, '{}_{}_BWA.log'.format(name, lane))
    nlog = os.path.join(errdir, '{}_{}_SAMTOOLS_SNAME.log'.format(name, lane))
    flog = os.path.join(errdir, '{}_{}_SAMTOOLS_FIX.log'.format(name, lane))
    rlog = os.path.join(errdir, '{}_{}_SAMTOOLS_ADDREPLACERG.log'.format(name, lane))
    olog = os.path.join(errdir, '{}_{}_SAMTOOLS_SORT.log'.format(name, lane))
    mlog = os.path.join(errdir, '{}_{}_SAMTOOLS_MARK.log'.format(name, lane))
    with open(blog, 'w') as ebwa, open(nlog, 'w') as enam, open(flog, 'w') as efix, \
         open(rlog, 'w') as earg, open(olog, 'w') as osor, open(mlog, 'w') as emar:
            pbwa = sp.Popen(shlex.split(curr_cmd_bwa), stdout=sp.PIPE, stderr=ebwa)
            pnam = sp.Popen(shlex.split(curr_cmd_nam), stdin=pbwa.stdout, stdout=sp.PIPE, stderr=enam)
            pfix = sp.Popen(shlex.split(curr_cmd_fix), stdin=pnam.stdout, stdout=sp.PIPE, stderr=efix)
            parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pfix.stdout, stdout=sp.PIPE, stderr=earg)
            psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=osor)
            pmar = sp.Popen(shlex.split(curr_cmd_mar), stdin=psor.stdout, stdout=sp.PIPE, stderr=emar)
            rcodes = map(lambda p : p.wait(), [pbwa, pnam, pfix, parg, psor, pmar])
    return '{}_{}'.format(name, lane), bam, check_rcodes(rcodes, [blog, nlog, flog, rlog, olog, mlog], [nam_tmp, sor_tmp, mar_tmp])


def barcode(files, names, barcodes, lanes, tmpdir, errdir, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_barcoding, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c is None else error_code(e, c, errdir))
    bams = [b for e, b, c in pool.imap_unordered(barcoding, jobs) if progress(e, c)]
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
    curr_cmd_arg = cmd_arg.format(fil, 'CB:Z:{}-{}'.format(barcode, lane), '{}'.format(name), '{}'.format(lane), bam)
    rlog = os.path.join(errdir, '{}_{}_SAMTOOLS_ADDREPLACERG.log'.format(name, lane))
    with open(rlog, 'w') as earg:
        parg = sp.Popen(shlex.split(curr_cmd_arg), stdout=sp.PIPE, stderr=earg)
        rcodes = map(lambda p : p.wait(), [parg])
    return '{}_{}'.format(name, lane), bam, check_rcodes(rcodes, [rlog])


def barcode_marked(files, names, barcodes, lanes, tmpdir, errdir, samtools, J):
    jobs = zip(files, names, barcodes, lanes)
    bar = ProgressBar(total=len(files), length=30, verbose=False)
    initargs = (tmpdir, errdir, samtools)
    pool = mp.Pool(processes=min(J, len(names)), initializer=init_barcoding_marked, initargs=initargs)
    progress = (lambda e, c : bar.progress(advance=True, msg="{}".format(e)) if c is None else error_code(e, c, errdir))
    bams = [b for e, b, c in pool.imap_unordered(barcoding_marked, jobs) if progress(e, c)]
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
    nlog = os.path.join(errdir, '{}_{}_SAMTOOLS_SNAME.log'.format(name, lane))
    flog = os.path.join(errdir, '{}_{}_SAMTOOLS_FIX.log'.format(name, lane))
    rlog = os.path.join(errdir, '{}_{}_SAMTOOLS_ADDREPLACERG.log'.format(name, lane))
    olog = os.path.join(errdir, '{}_{}_SAMTOOLS_SORT.log'.format(name, lane))
    mlog = os.path.join(errdir, '{}_{}_SAMTOOLS_MARK.log'.format(name, lane))
    with open(nlog, 'w') as enam, open(flog, 'w') as efix, open(rlog, 'w') as earg, open(olog, 'w') as osor, open(mlog, 'w') as emar:
        pnam = sp.Popen(shlex.split(curr_cmd_nam), stdout=sp.PIPE, stderr=enam)
        pfix = sp.Popen(shlex.split(curr_cmd_fix), stdin=pnam.stdout, stdout=sp.PIPE, stderr=efix)
        parg = sp.Popen(shlex.split(curr_cmd_arg), stdin=pfix.stdout, stdout=sp.PIPE, stderr=earg)
        psor = sp.Popen(shlex.split(curr_cmd_sor), stdin=parg.stdout, stdout=sp.PIPE, stderr=osor)    
        pmar = sp.Popen(shlex.split(curr_cmd_mar), stdin=psor.stdout, stdout=sp.PIPE, stderr=emar)
        rcodes = map(lambda p : p.wait(), [pnam, pfix, parg, psor, pmar])
    return '{}_{}'.format(name, lane), bam, check_rcodes(rcodes, [nlog, flog, rlog, olog, mlog], [nam_tmp, sor_tmp, mar_tmp])


def merging(samtools, bams, jobs, tmpdir, errdir):
    barcoded = os.path.join(tmpdir, 'barcodedcells.bam')
    cellbams = os.path.join(tmpdir, 'cellbams.tsv')
    with open(cellbams, 'w') as o:
        o.write('\n'.join(bams))
    cmd = '{} merge {} -b {} -@ {}'.format(samtools, barcoded, cellbams, jobs)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise ValueError(error('Merging failed with messages:\n{}\n{}\n'.format(stdout, stderr)))
    return barcoded


def indexing(samtools, jobs, tmpdir, errdir, barcoded, output):
    shutil.move(barcoded, output)
    cmd = '{} index {} -@ {}'.format(samtools, output, jobs)
    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise ValueError(error('Indexing failed with messages:\n{}\n{}\n'.format(stdout, stderr)))
    return


def error_code(e, c, errdir):
    raise ValueError(error('Commands failed for {} with error:\n{}\nCheck the errors in the corresponding log files in:\n{}'.format(e, c, errdir)))


def check_rcodes(rcodes, logs, tdir, cmd=None):
    res = ''
    for code, log in zip(rcodes, logs):
        if code != 0:
            if cmd is not None:
                return cmd
            with open(log, 'r') as i:
                res += 'From {}:\n'.format(log) + '\n'.join(i.readlines()) + '\n'  
    if res == '':
        map(lambda f : os.remove(f), logs)
        map(lambda d : shutil.rmtree(d), tdir)
        return None
    else:
        return res


if __name__ == '__main__':
    main()
