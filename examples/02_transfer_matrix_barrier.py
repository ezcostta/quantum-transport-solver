import numpy as np

from qtransport.regions import OneDimensionalRegion
from qtransport.transfer import OneDimensionalTransferMatrixSolver


def main():
    region = OneDimensionalRegion(
        onsite=[0.5] * 20,
        hopping=1.0,
    )

    solver = OneDimensionalTransferMatrixSolver(region)

    energies = np.linspace(-1.9, 1.9, 101)

    print("# energy transmission")
    for energy in energies:
        result = solver.solve(float(energy))
        print(f"{energy:.8f} {result.transmission:.12f}")


if __name__ == "__main__":
    main()
