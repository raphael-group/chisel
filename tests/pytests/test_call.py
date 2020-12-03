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
               'f5b1005a730048efd8213465ed0b8398'
