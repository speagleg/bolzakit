import hashlib

import pytest


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="run slow tests")
    parser.addoption("--runnetwork", action="store_true", default=False,
                     help="run network tests even if the Strohmaier-Uski file is not cached (may download it)")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow test (needs --runslow)")
    config.addinivalue_line("markers", "network: needs --runnetwork, or the Strohmaier-Uski file already cached "
                                       "with the right sha256")


def su_cached() -> bool:
    """True iff the Strohmaier-Uski file is in the cache and passes its sha256 check. Never downloads."""
    from bolzakit.spectral import eigen as E
    p = E.cache_dir() / E.SU_NAME
    return p.exists() and hashlib.sha256(p.read_bytes()).hexdigest() == E.SU_SHA256


def pytest_collection_modifyitems(config, items):
    skip_slow = pytest.mark.skip(reason="needs --runslow")
    skip_net = pytest.mark.skip(reason="needs --runnetwork (Strohmaier-Uski file not cached)")
    run_net = config.getoption("--runnetwork") or su_cached()
    for it in items:
        if "slow" in it.keywords and not config.getoption("--runslow"):
            it.add_marker(skip_slow)
        if "network" in it.keywords and not run_net:
            it.add_marker(skip_net)
