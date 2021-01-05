import os
import hashlib
import shutil
import tempfile
from chisel.Cloner import main

this_dir = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(this_dir, 'data')


def test_clone():
    with tempfile.NamedTemporaryFile('w') as f:
        # The Cloner.main function overwrites the input file!
        # Make a temporary copy for the purpose of testing the call.
        with tempfile.NamedTemporaryFile('w') as f_calls:
            input_file = os.path.join(DATA_FOLDER, 'output', 'calls.tsv')

            shutil.copy(input_file, f_calls.name)
            main(args=[
                f_calls.name,
                '--seed', '12'
            ], stdout_file=f.name)

            assert hashlib.md5(open(f.name, 'rb').read()).hexdigest() == \
                   '2b94dde27bc7f5e2c9215dc1f890f5f7'
