import numpy as np
from qtransport.leads import OneDimensionalLead


class OneDimensionalScatteringSolver:
    """
    Scattering solver for a finite 1D tight-binding chain
    connected to identical left/right leads.
    """

    def __init__(self, onsite_region, hopping: float = 1.0):
        self.onsite_region = np.array(onsite_region, dtype=float)
        self.hopping = hopping
        self.lead = OneDimensionalLead(hopping=hopping, onsite=0.0)

    def solve(self, energy: float):
        """
        Solve for reflection and transmission amplitudes.

        Unknowns:
            psi_0, ..., psi_{N-1}, r, tau

        Incoming wave from left:
            psi_j = exp(ikj) + r exp(-ikj), j <= -1

        Transmitted wave to right:
            psi_j = tau exp(ikj), j >= N
        """

        N = len(self.onsite_region)
        k = self.lead.wave_number(energy)
        t = self.hopping

        number_unknowns = N + 2
        A = np.zeros((N + 2, number_unknowns), dtype=complex)
        b = np.zeros(N + 2, dtype=complex)

        r_index = N
        tau_index = N + 1

        # Equation at left boundary site j = -1
        # E psi_-1 = -t psi_-2 - t psi_0
        psi_m1_in = np.exp(-1j * k)
        psi_m2_in = np.exp(-2j * k)

        A[0, 0] = -t
        A[0, r_index] = energy * np.exp(1j * k) + t * np.exp(2j * k)
        b[0] = -(energy * psi_m1_in + t * psi_m2_in)

        # Equations inside the scattering region: j = 0,...,N-1
        for j in range(N):
            row = j + 1
            A[row, j] = energy - self.onsite_region[j]

            if j > 0:
                A[row, j - 1] = t
            else:
                # Coupling to psi_-1 = incoming + reflected
                A[row, r_index] = t * np.exp(1j * k)
                b[row] = -t * np.exp(-1j * k)

            if j < N - 1:
                A[row, j + 1] = t
            else:
                # Coupling to psi_N = tau exp(ikN)
                A[row, tau_index] = t * np.exp(1j * k * N)

        # Equation at right boundary site j = N
        # E psi_N = -t psi_{N-1} - t psi_{N+1}
        row = N + 1
        A[row, N - 1] = -t
        A[row, tau_index] = energy * np.exp(1j * k * N) + t * np.exp(1j * k * (N + 1))

        x = np.linalg.solve(A, b)

        psi_region = x[:N]
        r = x[r_index]
        tau = x[tau_index]

        return {
            "energy": energy,
            "k": k,
            "psi_region": psi_region,
            "r": r,
            "t": tau,
        }
