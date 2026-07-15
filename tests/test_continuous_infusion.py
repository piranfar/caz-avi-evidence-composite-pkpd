import numpy as np
import pytest

from cazavi_model.continuous_infusion import (
    steady_state_concentration,
    unbound_concentration,
)


def test_steady_state_concentration_scalar() -> None:
    result = steady_state_concentration(500.0, 10.0)
    assert result == pytest.approx(50.0)


def test_steady_state_concentration_vectorized() -> None:
    result = steady_state_concentration(
        np.array([500.0, 1000.0]),
        np.array([10.0, 20.0]),
    )
    np.testing.assert_allclose(result, np.array([50.0, 50.0]))


def test_invalid_clearance_rejected() -> None:
    with pytest.raises(ValueError, match="strictly positive"):
        steady_state_concentration(500.0, 0.0)


def test_unbound_concentration() -> None:
    result = unbound_concentration(50.0, 0.9)
    assert result == pytest.approx(45.0)


def test_invalid_free_fraction_rejected() -> None:
    with pytest.raises(ValueError, match="between 0 and 1"):
        unbound_concentration(50.0, 1.1)
