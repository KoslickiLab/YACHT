import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="Run slow tests"
    )

def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runslow"):
        # Skip slow tests unless --runslow is specified
        skip_slow = pytest.mark.skip(reason="Need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)