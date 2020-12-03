import os
import pytest


@pytest.fixture(scope='module')
def input_folder():
    return os.getenv('TEST_DIRECTORY', os.path.join(os.path.dirname(__file__), 'data', 'input'))
