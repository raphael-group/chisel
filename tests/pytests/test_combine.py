import os
import hashlib
import tempfile
from chisel.Combiner import main

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_combine():
    with tempfile.NamedTemporaryFile('w') as f:
        main(args=[
            '-r', os.path.join(DATA_FOLDER, 'output', 'rdr.tsv'),
            '-b', os.path.join(DATA_FOLDER, 'output', 'baf.tsv'),
            '-j', '1',
            '-s', '12'
        ], stdout_file=f.name)

        assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
               'e6860842b280054321588d2153b50b94'
