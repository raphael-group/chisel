import os
import hashlib
import tempfile
from chisel.Caller import main

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_call():
    with tempfile.NamedTemporaryFile('w') as f:
        main(args=[
            os.path.join(DATA_FOLDER, 'output', 'combo.tsv'),
            '-j', '1',
            '--seed', '12'
        ], stdout_file=f.name)

        assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
               'be0b7472a1c8d6f8f96b9aa0fb8df3c9'
