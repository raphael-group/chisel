import os
import hashlib
import tempfile
from chisel.RDREstimator import main


def test_rdr(input_folder):
    with tempfile.NamedTemporaryFile('w') as f:
        main(args=[
            '-n', os.path.join(input_folder, 'normal.bam'),
            '-t', os.path.join(input_folder, 'cells.bam'),
            '-r', os.path.join(input_folder, 'hg19.fa'),
            '-b', '5Mb',
            '-m', '100000',
            '-c', 'chr6',
            '-j', '1'
        ], stdout_file=f.name)

        assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
               'eaf470105df0dcd8be7a995f1d6a8525'
