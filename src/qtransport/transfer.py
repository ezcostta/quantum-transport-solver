import numpy as np

from qtransport.leads import OneDimensionalLead
from qtransport.regions import OneDimensionalRegion
from qtransport.results import ScatteringResult


class OneDimensionalTransferMatrixSolver:
    def __init__(
        self,
        region: OneDimensionalRegion,
        left_lead: OneDimensionalLead | None = None,
        right_lead: OneDimensionalLead | None = None,
    ):
        self.region = region
        self.left_lead = left_lead or OneDimensionalLead(
            hopping=region.hopping,
            onsite=0.0,
        )
        self.right_lead = right_lead or OneDimensionalLead(
            hopping=region.hopping,
            onsite=0.0,
        )

        if self.left_lead.hopping != self.right_lead.hopping:
            raise NotImplementedError("Different lead hoppings are not supported yet.")

        if self.left_lead.onsite != self.right_lead.onsite:
            raise NotImplementedError(
                "Different lead onsite energies are not supported yet."
            )

    def local_matrix(self, energy: float, onsite: float) -> np.ndarray:
        t = self.region.hopping

        return np.array(
            [
                [-(energy - onsite) / t, -1.0],
                [1.0, 0.0],
            ],
            dtype=complex,
        )

    def total_matrix(self, energy: float) -> np.ndarray:
        total = np.eye(2, dtype=complex)

        for onsite in self.region.onsite:
            total = self.local_matrix(energy, onsite) @ total

        return total

    def solve(self, energy: float) -> ScatteringResult:
        k = self.left_lead.wave_number(energy)
        M = self.total_matrix(energy)

        N = self.region.size

        # Left side:
        # [psi_0, psi_-1]^T =
        # [1 + r, exp(-ik) + r exp(ik)]^T
        #
        # Right side:
        # [psi_N, psi_{N-1}]^T =
        # tau [exp(ikN), exp(ik(N-1))]^T
        #
        # Transfer relation:
        # [psi_N, psi_{N-1}]^T = M [psi_0, psi_-1]^T

        left_in = np.array([1.0, np.exp(-1j * k)], dtype=complex)
        left_ref = np.array([1.0, np.exp(1j * k)], dtype=complex)
        right_out = np.array(
            [np.exp(1j * k * N), np.exp(1j * k * (N - 1))],
            dtype=complex,
        )

        # M(left_in + r left_ref) = tau right_out
        # r M left_ref - tau right_out = - M left_in

        A = np.column_stack((M @ left_ref, -right_out))
        b = -(M @ left_in)

        r, tau = np.linalg.solve(A, b)

        return ScatteringResult(
            energy=energy,
            reflection_amplitude=r,
            transmission_amplitude=tau,
        )
