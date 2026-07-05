from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from qtransport.regions import OneDimensionalRegion
from qtransport.transfer import OneDimensionalTransferMatrixSolver


def main():
    region = OneDimensionalRegion(
        onsite=[0.5] * 20,
        hopping=1.0,
    )

    solver = OneDimensionalTransferMatrixSolver(region)

    energies = np.linspace(-1.9, 1.9, 500)
    transmissions = []

    for energy in energies:
        result = solver.solve(float(energy))
        transmissions.append(result.transmission)

    output_dir = Path("figures/generated")
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 4))
    plt.plot(energies, transmissions, linewidth=2)
    plt.xlabel(r"Energy $E/t$")
    plt.ylabel(r"Transmission $T(E)$")
    plt.title("Transmission through a finite 1D barrier")
    plt.ylim(-0.05, 1.05)
    plt.tight_layout()

    plt.savefig(output_dir / "02_transfer_matrix_barrier.pdf")
    plt.savefig(output_dir / "02_transfer_matrix_barrier.png", dpi=300)

    print(f"Saved figure to {output_dir}")


if __name__ == "__main__":
    main()
