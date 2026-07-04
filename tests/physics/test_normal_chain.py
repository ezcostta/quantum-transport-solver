import pytest

from qtransport.solver import OneDimensionalScatteringSolver
from qtransport.observables import (
    reflection_probability,
    transmission_probability,
)


def test_probability_conservation():

    solver = OneDimensionalScatteringSolver(
        onsite_region=[0.5, 0.5, 0.5],
        hopping=1.0,
    )

    result = solver.solve(0.3)

    R = reflection_probability(result)
    T = transmission_probability(result)

    assert abs(R + T - 1) < 1e-12


def test_uniform_chain():

    solver = OneDimensionalScatteringSolver(
        onsite_region=[0, 0, 0],
        hopping=1.0,
    )

    result = solver.solve(0.3)

    R = reflection_probability(result)
    T = transmission_probability(result)

    assert R < 1e-12
    assert abs(T - 1) < 1e-12


def test_barrier():

    solver = OneDimensionalScatteringSolver(
        onsite_region=[2],
        hopping=1.0,
    )

    result = solver.solve(0.3)

    T = transmission_probability(result)

    assert 0 < T < 1


def test_energy_outside_band():

    solver = OneDimensionalScatteringSolver(
        onsite_region=[0],
        hopping=1.0,
    )

    with pytest.raises(ValueError):
        solver.solve(3.0)
