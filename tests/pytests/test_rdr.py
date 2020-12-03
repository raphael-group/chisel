import os
import hashlib
import tempfile
from chisel.RDREstimator import main

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_rdr():
    with tempfile.NamedTemporaryFile('w') as f:
        main(args=[
            '-n', os.path.join(DATA_FOLDER, 'normal.bam'),
            '-t', os.path.join(DATA_FOLDER, 'cells.bam'),
            '-r', os.path.join(DATA_FOLDER, 'hg19.fa'),
            '-b', '5Mb',
            '-m', '100000',
            '-c', 'chr6',
            '-j', '1',
            '-s', '/home/vineetb/.conda/envs/chisel/bin/samtools'
        ], stdout_file=f.name)

        assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
               'eaf470105df0dcd8be7a995f1d6a8525'
