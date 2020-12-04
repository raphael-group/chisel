import os
import hashlib
import tempfile
from chisel.BAFEstimator import main

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_baf(input_folder):
    with tempfile.NamedTemporaryFile('w') as f:
        # The complete chr6 phase data takes quite a while to process.
        # Here we test with just the first 1k rows of the phase data.
        with tempfile.NamedTemporaryFile('w') as phases_f:
            with open(os.path.join(input_folder, 'phases.tsv'), 'r') as phases_full:
                lines = [next(phases_full) for _ in range(1000)]
            phases_f.writelines(lines)
            phases_f.flush()
            main(args=[
                '-n', os.path.join(input_folder, 'normal.bam'),
                '-t', os.path.join(input_folder, 'cells.bam'),
                '-r', os.path.join(input_folder, 'hg19.fa'),
                '-j', '1',
                '-c', os.path.join(DATA_FOLDER, 'output', 'total.tsv'),
                '-l', phases_f.name
            ], stdout_file=f.name)

            assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
                   'f695990677d8488e82198d1b3c270a33'
