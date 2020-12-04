import os


def test_env():
    assert os.getenv('TEST_DIRECTORY') == '/home/runner/work/chisel/chisel/testdata'
