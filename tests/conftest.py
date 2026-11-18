import os, sys
import numpy as np
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cislab.constants import MU_EARTH_MOON


@pytest.fixture(scope="session")
def mu():
    return MU_EARTH_MOON


@pytest.fixture(scope="session")
def family(mu):
    """The Earth-Moon L2 northern halo family (built once per test session)."""
    from cislab.cr3bp.families import build_family
    return build_family(mu)
