#!/usr/bin/env python2.7

import os, sys
os.environ["OMP_NUM_THREADS"] = "1" 
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
import shlex
import argparse
import subprocess as sp

from multiprocessing import Lock, Value, Pool
from collections import Counter

import numpy as np
from scipy.stats import beta
from scipy.stats import binom
from scipy.special import binom as binom_coeff
from scipy.special import beta as beta_sp

from Utils import *


def parse_args(args):
    description = "Compute RDR from barcoded single-cell sequencing data."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-t","--tumor", required=True, type=str, help="Barcoded BAM file for tumor sample")
    parser.add_argument("-n","--normal", required=False, default=None, type=str, help="BAM file for matched normal sample")
    parser.add_argument("-l","--listphased", type=str, required=True, help="List of TSV phased genomic positions (TAB-seprated format '#CHR POS PHASE' where PHASE is either 0|1 or 1|0)")
    parser.add_argument("-r","--reference", type=str, required=True, help="Reference genome")
    parser.add_argument("-b","--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("-s","--samtools", required=False, default=None, type=str, help="Path to the directory to \"samtools\" executable, required in default mode (default: samtools is directly called as it is in user $PATH)")
    parser.add_argument("-j","--jobs", required=False, type=int, default=0, help="Number of parallele jobs to use (default: equal to number of available processors)")
    parser.add_argument("-c","--listcells", type=str, required=False, default=None, help="File where first column contains all the cells to consider (default: not used)")
    parser.add_argument("--cellprefix", type=str, required=False, default='CB:Z:', help="Prefix of cell barcode field in SAM format (default: CB:Z:)")
    parser.add_argument("--cellsuffix", type=str, required=False, default='', help="Suffix of cell barcode field in SAM format (default: none)")
    parser.add_argument("--rundir", type=str, required=False, default='./', help="Running directory (default: ./)")
    args = parser.parse_args(args)

    if not os.path.isfile(args.tumor):
        raise ValueError("Specified tumor does not exist!")
    if args.normal is not None and not os.path.isfile(args.normal):
        raise ValueError("Specified normal does not exist!")
    elif args.normal is None:
        log('A matched-normal sample has not been provided and the presence of the provided heterozygous germline SNPs will not be correspondigly assessed.', level='WARN')
    if not os.path.isfile(args.listphased):
        raise ValueError("List of phased postions does not exist!")
    if not os.path.isfile(args.reference):
        raise ValueError("Reference genome does not exist!")
    if args.listcells is not None and not os.path.isfile(args.listcells):
        raise ValueError("Specified list of cells does not exist!")
    if not os.path.isdir(args.rundir):
        raise ValueError("Running directory does not exists: {}".format(args.rundir))

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
        'tumor' : args.tumor,
        'normal' : args.normal,
        'phased' : args.listphased,
        'ref' : args.reference,
        'bcftools' : bcftools,
        'samtools' : samtools,
        'J' : args.jobs,
        'gamma' : 0.01,
        'seqerr' : 0.001,
        'list' : args.listcells,
        'cellprefix' : args.cellprefix,
        'cellsuffix' : args.cellsuffix,
        'rundir' : os.path.abspath(args.rundir)
    }


def main(args=None, stdout_file=None):
    log('Parsing and checking arguments')
    args = parse_args(args)
    log('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]), level='INFO')

    log('Extracting chromosomes')
    phased = read_phase(args['phased'])
    log('Chromosomes analyzed: {}'.format(','.join(sorted(phased, key=orderchrs))), level='INFO')
    log('Total number of given phased positions: {}'.format(sum(len(phased[c]) for c in phased)), level='INFO')

    if args['normal'] is not None:
        log('Counting phased SNPs in matched normal')
    else:
        log('Selecting phased SNPs')
    snps = selecting(args, phased)
    log('Number of selected heterozygous SNPs: {}'.format(sum(len(snps[c]) for c in snps)), level='INFO')
    # selfile = os.path.join(args['rundir'], 'selectedSNPs.tsv')
    # log('Writing selected phased SNPs in {}'.format(selfile), level='INFO')
    # with open(selfile, 'w') as o:
    #     det = (lambda c, o : [snps[c][o]['REFALT'][0], snps[c][o]['REFALT'][1], snps[c][o]['NORM'][0], snps[c][o]['NORM'][1], snps[c][o]['PHASE']])
    #     rec = (lambda c, o : '\t'.join(map(str, [c, o] + det(c, o))))
    #     o.write('\n'.join([rec(c, o) for c in sorted(snps, key=orderchrs) for o in sorted(snps[c])]) + '\n')

    log('Extracting SNP counts for all cells')
    abc = extracting(args, snps)

    if args['list']:
        log('Reading cell list')
        with open(args['list'], 'r') as i:
            cells = set(l.strip().split()[0].replace('-1', '') for l in i if len(l) > 1 and l[0] != '#')

    log('Writing A/B counts for selected phased SNPs across selected cells')
    form = (lambda c, o, s, REF, ALT : '\t'.join(map(str, [c, o, s] + ([REF, ALT] if snps[c][o]['PHASE'] == '0|1' else [ALT, REF]))))
    getREF = (lambda c, o, s : abc[c][o][s][snps[c][o]['REFALT'][0]])
    getALT = (lambda c, o, s : abc[c][o][s][snps[c][o]['REFALT'][1]])
    check = (lambda c, o, s : (getREF(c, o, s) + getALT(c, o, s)) >= 0)
    if args['list']:
        rec = (lambda c, o : (form(c, o, s, getREF(c, o, s), getALT(c, o, s)) for s in sorted(abc[c][o]) if check(c, o, s) and s in cells))
    else:
        rec = (lambda c, o : (form(c, o, s, getREF(c, o, s), getALT(c, o, s)) for s in sorted(abc[c][o]) if check(c, o, s)))
    #print '\n'.join([r for c in sorted(snps, key=orderchrs) for o in sorted(snps[c]) for r in rec(c, o)])

    if stdout_file is not None:
        with open(stdout_file, 'w') as stdout_f:
            for p in (r for c in sorted(snps, key=orderchrs) for o in sorted(snps[c]) for r in rec(c, o)):
                stdout_f.write("{}\n".format(p))
    else:
        for p in (r for c in sorted(snps, key=orderchrs) for o in sorted(snps[c]) for r in rec(c, o)):
            print p
    
    log('KTHXBYE')


def read_phase(f):
    phased = defaultdict(lambda : dict())
    with open(f, 'r') as i:
        for l in i:
            p = l.strip().split()
            if len(l) > 1 and p[0][0] != '#':
                zeroone = '0|1' in l
                onezero = '1|0' in l
                if zeroone or onezero:
                    if zeroone and onezero:
                        raise ValueError('Found a record in phased positions which contains both phases 0|1 and 1|0!')
                    if p[0] in phased[p[0]]:
                        raise ValueError('Found a duplicate phased position!')
                    phased[p[0]][int(p[1])] = '0|1' if zeroone else '1|0'
    return phased


def selecting(args, phased):
    jobs = [(c, os.path.join(args['rundir'], 'hetSNPs-{}.tmp'.format(c))) for c in phased]
    for c, f in jobs:
        with open(f, 'w') as o:
            o.write('\n'.join(['{}\t{}'.format(c, o) for o in sorted(phased[c])]) + '\n')

    lock = Lock()
    counter = Value('i', 0)
    bam = args['normal'] if args['normal'] is not None else args['tumor']
    initargs = (lock, counter, len(jobs), bam, args['bcftools'], args['ref'], args['gamma'], args['seqerr'], args['normal'] is not None)
    pool = Pool(processes=min(args['J'], len(phased)), initializer=init_selecting, initargs=initargs)

    uniq = (lambda l : reduce(lambda x, y : inupdate(x, y), [{e[0] : (e[1], e[2], e[3], e[4])} for e in reversed(l)]))
    refalt = {c : uniq(l) for c, l in pool.imap_unordered(counting_germinal, jobs) if len(l) > 0}
    pool.close()
    pool.join()
    form = (lambda c, o : {'REFALT' : tuple(refalt[c][o][:2]), 'PHASE' : phased[c][o], 'NORM' : tuple(refalt[c][o][2:])})
    return {c : {o : form(c, o) for o in refalt[c]} for c in refalt}


def init_selecting(lock, counter, l, bam, bcf, ref, _gamma, _seqerr, _isnormal):
    global bar, cmd_mpi, cmd_que, gamma, seqerr, isnormal
    bar = ProgressBar(total=l, length=40, lock=lock, counter=counter, verbose=False)
    cmd_mpi = '{} mpileup {} -R {} -f {} --skip-indels -a INFO/AD -Ou'.format(bcf, bam, "{}", ref)
    cmd_que = "{} query {}".format(bcf, "-f '%CHROM\t%POS\t%REF\t%ALT{0}\t%AD{0}\t%AD{1}\n' -i 'SUM(AD)<=10000'")
    gamma = _gamma
    seqerr = _seqerr
    isnormal = _isnormal


def counting_germinal(j):
    c, f = j
    bar.progress(advance=False, msg="Evaluating germinal SNPs on {}".format(c))
    mpi = sp.Popen(shlex.split(cmd_mpi.format(f)), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(shlex.split(cmd_que), stdin=mpi.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    if stderr:
        raise ValueError("ERROR: {}".format(stderr))
        sys.exit(1)
    os.remove(f)

    ACGT = {'A', 'C', 'G', 'T'}
    def asser(C):
        assert C == c
        return True
    if isnormal:
        check = (lambda C, O, ref, alt, A, B : ref in ACGT and alt in ACGT and isHet(int(A), int(B), gamma) and asser(C))
    else:
        check = (lambda C, O, ref, alt, A, B : ref in ACGT and alt in ACGT and isVar(int(A), int(B), gamma, seqerr) and asser(C))
    form = (lambda C, O, ref, alt, A, B : (int(O), ref, alt, A, B))

    bar.progress(advance=True, msg="Finalizing germinal SNPs on {}".format(c))
    return (c, [form(*l.split()) for l in stdout.strip().split('\n') if l != '' and check(*l.split())])


def isHet(countA, countB, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= 0.5 <= c_upper


def isVar(countA, countB, gamma, eps, lam=0.5):
    p0 = probX(False, min(countA, countB), countA + countB, lam, eps) # Absence variant
    p1 = probX(True, min(countA, countB), countA + countB, lam, eps) # Presence variant
    if p1 == float('inf'):
        return 0, 1
    return (p1 / (p0 + p1)) > 1.0 - gamma


def probX(X, V, T, lam, eps):
    if X:
        b = binom_coeff(T, V)
        if b == float('inf'): return float('inf')
        return lam * binom_coeff(T, V) * beta_sp(1 + V, 1 + T - V)
    else:
        return (1 - lam) * binom.pmf(V, T, eps)


def extracting(args, snps):
    countawk = os.path.join(args['rundir'], 'count.awk') #os.path.join(os.path.dirname(os.path.abspath(__file__)), 'count.awk')
    mkcount(countawk)
    ACGT = (lambda : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0, 'N' : 0})
    abc = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : ACGT())))
    abc = {c : {o : defaultdict(lambda : ACGT()) for o in snps[c]} for c in snps}
    bar = ProgressBar(total=sum(len(snps[chro]) for chro in snps), length=40, verbose=False)
            
    for chro in sorted(snps, key=orderchrs):
        jobs = [(chro, o) for o in snps[chro]]
        initargs = (args['tumor'], args['samtools'], countawk, args['cellprefix'], args['cellsuffix'])
        pool = Pool(processes=min(args['J'], len(jobs)), initializer=init_extracting, initargs=initargs)
        
        for c, o, l in pool.imap_unordered(counting_cell, jobs):
            if l != '':
                for a in l.strip().split('\n'):
                    e, al, count = tuple(a.split())
                    abc[c][o][e][al] += int(count)
            bar.progress(advance=True, msg="Extracted SNP {}:{}".format(c, o))
            
        pool.close()
        pool.join()

    os.remove(countawk)
    return {c : {o : dict(filter(lambda (e, al) : sum(al.values()) > 0, abc[c][o].items())) for o in abc[c]} for c in abc}


def init_extracting(_tumor, _sam, countawk, cellprefix, cellsuffix):
    global cmd_sam, cmd_awk
    cmd_sam = "{} view -F 1796 -q 13 {} {}:{}-{}".format(_sam, _tumor, '{}', '{}', '{}')
    cmd_awk = 'awk -v TAG="{}" -v CELLPREFIX="{}" -v CELLSUFFIX="{}" -f {}'.format('{}', cellprefix, cellsuffix, countawk)


def counting_cell(job):
    sam = sp.Popen(shlex.split(cmd_sam.format(job[0], job[1], job[1])), stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = sp.Popen(shlex.split(cmd_awk.format(job[1])), stdin=sam.stdout, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    return (job[0], job[1], stdout)


def mkcount(f):
    with open(f, 'w') as o:
        o.write('#!/usr/bin/awk\n\n')
        o.write('BEGIN{}\n')
        o.write('{\n')
        o.write('    if ( match($0, CELLPREFIX"[ACGT]+"CELLSUFFIX) )\n')
        o.write('    {\n')
        o.write('        REF = $4 - 1;\n')
        o.write('        QUE = 0;\n')
        o.write('        CIG = $6;\n')
        o.write('        CEL = substr($0, RSTART+5, RLENGTH-5);\n')
        o.write('        while( match(CIG, /^[[:digit:]]+/) )\n')
        o.write('        {\n')
        o.write('            N = substr(CIG, RSTART, RLENGTH);\n')
        o.write('            CIG = substr(CIG, RSTART+RLENGTH);\n')
        o.write('            if( match(CIG, /^[MIDNSHP=X]/) )\n')
        o.write('            {\n')
        o.write('                C = substr(CIG, RSTART, RLENGTH);\n')
        o.write('                CIG = substr(CIG, RSTART+RLENGTH);\n')
        o.write('                if (C == "M" || C == "=" || C == "X")\n')
        o.write('                {\n')
        o.write('                    REF += N;\n')
        o.write('                    QUE += N;\n')
        o.write('                    if (TAG <= REF)\n')
        o.write('                    {\n')
        o.write('                        X[CEL, substr($10, QUE - REF + TAG, 1)]++;\n')
        o.write('                        next;\n')
        o.write('                    };\n')
        o.write('                } else if (C == "D" || C == "N")\n')
        o.write('                {\n')
        o.write('                    REF += N;\n')
        o.write('                    if (TAG <= REF)\n')
        o.write('                    {\n')
        o.write('                        X[CEL, "N"]++;\n')
        o.write('                        next;\n')
        o.write('                    }\n')
        o.write('                } else if (C == "I" || C == "S")\n')
        o.write('                {\n')
        o.write('                    QUE += N;\n')
        o.write('                };\n')
        o.write('            };\n')
        o.write('        };\n')
        o.write('    };\n')
        o.write('}\n')
        o.write('END{ for (p in X) { split(p, x, SUBSEP); print x[1], x[2], X[x[1], x[2]] } }\n')
    return


if __name__ == '__main__':
    main()
