import pytest


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="run slow tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow test (needs --runslow)")
    config.addinivalue_line("markers", "network: needs network or the cached Strohmaier-Uski file")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="needs --runslow")
    for it in items:
        if "slow" in it.keywords:
            it.add_marker(skip)
