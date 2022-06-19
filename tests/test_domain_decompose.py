import pytest
from domain_decompose import DomainDecompose


@pytest.fixture
def test_domain_decompose():
    return DomainDecompose()
