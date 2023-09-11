import pytest

@pytest.fixture
def tf_optimizer():
    return 1

#this needs to be moved into test_installation if change to pytest
def test_tfoptimizer(tf_optimizer):
    print(tf_optimizer)
