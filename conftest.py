# 
# conftest.py                                                               
# 
# D. Clarke 
# 
# Configuration file for pytest 
# 

import pytest

@pytest.fixture(autouse=True)
def executeInPlace(request, monkeypatch):
    """ This executes each test in the directory that test sits in. """
    monkeypatch.chdir(request.fspath.dirname)

