import numpy as np

from qtransport.leads import OneDimensionalLead
from qtransport.regions import OneDimensionalRegion
from qtransport.results import ScatteringResult


class OneDimensionalTransferMatrixSolver:
    """
    Transfer-matrix solver for a one-dimensional nearest-neighbor
    tight-binding chain.

    The local propagation equation is

        Phi_{j+1} = M_j Phi_j,

    where

        Phi_j = (psi_j, psi_{j-1})^T.

    This implements the transfer-matrix method described in the
    Theory Manual, Chapter 2.
    """

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

        self._validate_leads()

    def _validate_leads(self) -> None:
        if self.left_lead.hopping != self.right_lead.hopping:
            raise NotImplementedError("Different lead hoppings are not supported yet.")

        if self.left_lead.onsite != self.right_lead.onsite:
            raise NotImplementedError(
                "Different lead onsite energies are not supported yet."
            )

    def local_matrix(self, energy: float, onsite: float) -> np.ndarray:
        """
        Construct the local transfer matrix M_j.

        The tight-binding equation

            (E - eps_j) psi_j + t psi_{j-1} + t psi_{j+1} = 0

        is rewritten as

            psi_{j+1} = -((E - eps_j)/t) psi_j - psi_{j-1}.

        Therefore,

            M_j = [[-(E - eps_j)/t, -1],
                   [1,                 0]].
        """
        hopping = self.region.hopping

        return np.array(
            [
                [-(energy - onsite) / hopping, -1.0],
                [1.0, 0.0],
            ],
            dtype=complex,
        )

    def total_matrix(self, energy: float) -> np.ndarray:
        """
        Construct the total transfer matrix of the scattering region.

        If the scattering region contains N sites, then

            M = M_{N-1} M_{N-2} ... M_1 M_0.

        The multiplication order matters because transfer matrices
        generally do not commute.
        """
        total = np.eye(2, dtype=complex)

        for onsite in self.region.onsite:
            total = self.local_matrix(energy, onsite) @ total

        return total

    def incoming_vector(self, wave_number: float) -> np.ndarray:
        """
        Return the incoming basis vector at the left boundary.

        Phi_in = (1, exp(-ik))^T.
        """
        return np.array(
            [1.0, np.exp(-1j * wave_number)],
            dtype=complex,
        )

    def reflected_vector(self, wave_number: float) -> np.ndarray:
        """
        Return the reflected basis vector at the left boundary.

        Phi_ref = (1, exp(ik))^T.
        """
        return np.array(
            [1.0, np.exp(1j * wave_number)],
            dtype=complex,
        )

    def transmitted_vector(self, wave_number: float) -> np.ndarray:
        """
        Return the outgoing transmitted basis vector at the right boundary.

        Phi_out = (exp(ikN), exp(ik(N-1)))^T.
        """
        number_of_sites = self.region.size

        return np.array(
            [
                np.exp(1j * wave_number * number_of_sites),
                np.exp(1j * wave_number * (number_of_sites - 1)),
            ],
            dtype=complex,
        )

    def solve(self, energy: float) -> ScatteringResult:
        """
        Solve the scattering problem at a fixed energy.

        The total transfer matrix satisfies

            Phi_N = M Phi_0.

        Using

            Phi_0 = Phi_in + r Phi_ref

        and

            Phi_N = tau Phi_out,

        we obtain

            r M Phi_ref - tau Phi_out = - M Phi_in.

        This is a 2x2 linear system for r and tau.
        """
        wave_number = self.left_lead.wave_number(energy)
        total_matrix = self.total_matrix(energy)

        incoming = self.incoming_vector(wave_number)
        reflected = self.reflected_vector(wave_number)
        transmitted = self.transmitted_vector(wave_number)

        coefficient_matrix = np.column_stack(
            (
                total_matrix @ reflected,
                -transmitted,
            )
        )

        right_hand_side = -(total_matrix @ incoming)

        reflection_amplitude, transmission_amplitude = np.linalg.solve(
            coefficient_matrix,
            right_hand_side,
        )

        return ScatteringResult(
            energy=energy,
            reflection_amplitude=reflection_amplitude,
            transmission_amplitude=transmission_amplitude,
        )
