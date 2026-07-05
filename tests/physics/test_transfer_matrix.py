import numpy as np

from qtransport.regions import OneDimensionalRegion
from qtransport.solver import OneDimensionalScatteringSolver
from qtransport.transfer import OneDimensionalTransferMatrixSolver
from qtransport.observables import transmission_probability


def test_transfer_matrix_uniform_chain():
    region = OneDimensionalRegion(onsite=[0.0, 0.0, 0.0], hopping=1.0)
    solver = OneDimensionalTransferMatrixSolver(region)

    result = solver.solve(energy=0.3)

    assert result.reflection < 1e-12
    assert abs(result.transmission - 1.0) < 1e-12


def test_transfer_matrix_probability_conservation():
    region = OneDimensionalRegion(onsite=[0.5, 0.5, 0.5], hopping=1.0)
    solver = OneDimensionalTransferMatrixSolver(region)

    result = solver.solve(energy=0.3)

    assert result.conservation_error < 1e-12


def test_transfer_matrix_matches_boundary_solver():
    onsite = [0.5, 0.5, 0.5, 0.5]
    energy = 0.3

    region = OneDimensionalRegion(onsite=onsite, hopping=1.0)
    tm_solver = OneDimensionalTransferMatrixSolver(region)

    bm_solver = OneDimensionalScatteringSolver(
        onsite_region=onsite,
        hopping=1.0,
    )

    tm_result = tm_solver.solve(energy)
    bm_result = bm_solver.solve(energy)

    assert np.isclose(
        tm_result.transmission,
        transmission_probability(bm_result),
        atol=1e-12,
    )
