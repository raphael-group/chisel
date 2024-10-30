#!/usr/bin/env python2.7

import os, sys
import sys
import datetime
import subprocess as sp
import multiprocessing as mp
import shlex
import datetime
import re

from collections import defaultdict
from contextlib import contextmanager


argmax = (lambda d : max(d, key=(lambda x : d[x])))


argmin = (lambda d : min(d, key=(lambda x : d[x])))


def log(m):
    sys.stderr.write('# ' + m + '\n')

checkchrs = (lambda digit_x, x : digit_x if len(digit_x) > 0 else 
                                (23 if 'X' in x else 
                                (24 if 'Y' in x else 
                                (25 if 'M' in x else 26))))
orderchrs = (lambda x : int(checkchrs(''.join([l for l in x if l.isdigit()]), x)))


def inupdate(a, b):
    a.update(b)
    return a


def indices_to_condensed(i, j, n):
    assert i != j, "ERROR: equal indices cannot be transformed into condensed indices!"
    if i < j:
        i, j = j, i
    return n*j - j*(j+1)/2 + i - 1 - j

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

@contextmanager
def stdout_redirected(to=os.devnull):
    fd = sys.stdout.fileno()

    def _redirect_stdout(to):
        sys.stdout.close()
        os.dup2(to.fileno(), fd)
        sys.stdout = os.fdopen(fd, 'w')

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'w') as file:
            _redirect_stdout(to=file)
        try:
            yield
        finally:
            _redirect_stdout(to=old_stdout)


class ProgressBar:

    def __init__(self, total, length, lock=None, counter=0, verbose=False, decimals=1, fill=unichr(9608), prefix = 'Progress:', suffix = 'Complete'):
        self.total = total
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.prefix = prefix
        self.suffix = suffix
        self.lock = lock
        self.counter = counter
        assert lock is not None or counter == 0
        self.verbose = verbose

    def progress(self, advance=True, msg=""):
        if self.lock is None:
            self.progress_unlocked(advance, msg)
        else:
            self.progress_locked(advance, msg)
        return True

    def progress_unlocked(self, advance, msg):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            self.counter += 1
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.counter / float(self.total)))
        filledLength = int(self.length * self.counter // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix)
        msg = '[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + msg
        if not self.verbose:
            toprint = rewind + result + " [%s]" % (msg)
        else:
            toprint = rewind + msg + "\n" + result
        write(toprint.encode('utf-8'))
        flush()
        if self.counter == self.total:
            write("\n")
            flush()

    def progress_locked(self, advance, msg):
        flush = sys.stderr.flush
        write = sys.stderr.write
        if advance:
            with self.counter.get_lock():
                self.counter.value += 1
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.counter.value / float(self.total)))
        filledLength = int(self.length * self.counter.value // self.total)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        rewind = '\x1b[2K\r'
        result = '%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix)
        msg = '[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + msg
        if not self.verbose:
            toprint = rewind + result + " [%s]" % (msg)
        else:
            toprint = rewind + msg + "\n" + result
        with self.lock:
            write(toprint.encode('utf-8'))
            flush()
            if self.counter.value == self.total:
                write("\n")
                flush()


def log(msg, level='STEP', lock=None):
    timestamp = '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
    if level == "STEP":
        if lock is None:
            sys.stderr.write("{}{}[{}]{}{}\n".format(bcolors.BOLD, bcolors.HEADER, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}{}[{}]{}{}\n".format(bcolors.BOLD, bcolors.HEADER, timestamp, msg, bcolors.ENDC))
    elif level == "INFO":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}\n".format(bcolors.OKGREEN, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}\n".format(bcolors.OKGREEN, timestamp, msg, bcolors.ENDC))
    elif level == "WARN":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}\n".format(bcolors.WARNING, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}\n".format(bcolors.WARNING, timestamp, msg, bcolors.ENDC))
    elif level == "PROGRESS":
        if lock is None:
            sys.stderr.write("{}{}[{}]{}{}\n".format(bcolors.UNDERLINE, bcolors.BBLUE, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}\n".format(bcolors.BBLUE, timestamp, msg, bcolors.ENDC))
    elif level == "ERROR":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}\n".format(bcolors.FAIL, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}\n".format(bcolors.FAIL, timestamp, msg, bcolors.ENDC))
    else:
        if lock is None:
            sys.stderr.write("{}\n".format(msg))
        else:
            with lock: sys.stderr.write("{}\n".format(msg))


def runcmd(cmd, xdir, out=None, log="log"):
    j = os.path.join
    tmp = log + '_TMP'
    sout = open(j(xdir, out) if out is not None else os.devnull, 'w')
    with open(j(xdir, tmp), 'w') as serr:
        proc = sp.Popen(shlex.split(cmd), stdout=sout, stderr=sp.PIPE)
        for line in iter(lambda : proc.stderr.read(1), ''):
            sys.stderr.write(line)
            serr.write(line)
    sout.flush()
    sout.close()

    with open(j(xdir, tmp), 'r') as i:
        with open(j(xdir, log), 'w') as o:
            for l in i:
                if 'Progress' not in l:
                    o.write(re.sub(r'\033\[[0-9]*m', '', l))
    os.remove(j(xdir, tmp))


def error(msg):
    log(msg=msg, level="ERROR")
    sys.exit(0)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    BBLUE = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
