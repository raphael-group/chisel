import os
import hashlib
import tempfile
from chisel.BAFEstimator import main

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_baf():
    with tempfile.NamedTemporaryFile('w') as f:
        main(args=[
            '-n', os.path.join(DATA_FOLDER, 'normal.bam'),
            '-t', os.path.join(DATA_FOLDER, 'cells.bam'),
            '-r', os.path.join(DATA_FOLDER, 'hg19.fa'),
            '-j', '1',
            '-c', os.path.join(DATA_FOLDER, 'output', 'total.tsv'),
            '-l', os.path.join(DATA_FOLDER, 'phasing.tsv'),
            '-s', '/home/vineetb/.conda/envs/chisel/bin/samtools',
            '-b', '/home/vineetb/.conda/envs/chisel/bin/bcftools'
        ], stdout_file=f.name)

        assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
               'f695990677d8488e82198d1b3c270a33'
