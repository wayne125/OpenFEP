from openfep.lambda_scheduler.presets import rbfe_lambdas
from openfep.constants import RBFE_DEFAULT_LAMBDAS
import pytest


def test_rbfe_lambdas_default_count():
    lams = rbfe_lambdas()
    assert len(lams) == RBFE_DEFAULT_LAMBDAS


def test_rbfe_lambdas_endpoints():
    lams = rbfe_lambdas()
    assert lams[0] == 0.0
    assert lams[-1] == 1.0


def test_rbfe_lambdas_custom_n():
    lams = rbfe_lambdas(n=5)
    assert len(lams) == 5
    assert lams[0] == 0.0
    assert lams[-1] == 1.0


def test_rbfe_lambdas_spacing():
    lams = rbfe_lambdas(n=3)
    assert abs(lams[1] - 0.5) < 1e-9


def test_rbfe_lambdas_too_small():
    with pytest.raises(ValueError):
        rbfe_lambdas(n=1)
