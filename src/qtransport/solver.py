import numpy as np

from qtransport.leads import OneDimensionalLead
from qtransport.results import ScatteringResult


class OneDimensionalScatteringSolver:
    """
    Boundary-matching solver for a finite one-dimensional tight-binding
    chain connected to identical semi-infinite leads.
    """

    def __init__(self, onsite_region, hopping: float = 1.0):
        self.onsite_region = np.array(onsite_region, dtype=float)
        self.hopping = hopping
        self.lead = OneDimensionalLead(hopping=hopping, onsite=0.0)

    def solve(self, energy: float) -> ScatteringResult:
        """
        Solve for reflection and transmission amplitudes using direct
        boundary matching.
        """

        number_of_sites = len(self.onsite_region)
        wave_number = self.lead.wave_number(energy)
        hopping = self.hopping

        number_unknowns = number_of_sites + 2
        matrix = np.zeros((number_of_sites + 2, number_unknowns), dtype=complex)
        rhs = np.zeros(number_of_sites + 2, dtype=complex)

        reflection_index = number_of_sites
        transmission_index = number_of_sites + 1

        psi_minus_one_in = np.exp(-1j * wave_number)
        psi_minus_two_in = np.exp(-2j * wave_number)

        # Left lead boundary equation at j = -1.
        matrix[0, 0] = hopping
        matrix[0, reflection_index] = energy * np.exp(
            1j * wave_number
        ) + hopping * np.exp(2j * wave_number)
        rhs[0] = -(energy * psi_minus_one_in + hopping * psi_minus_two_in)

        # Scattering region equations.
        for site in range(number_of_sites):
            row = site + 1

            matrix[row, site] = energy - self.onsite_region[site]

            if site > 0:
                matrix[row, site - 1] = hopping
            else:
                matrix[row, reflection_index] = hopping * np.exp(1j * wave_number)
                rhs[row] = -hopping * np.exp(-1j * wave_number)

            if site < number_of_sites - 1:
                matrix[row, site + 1] = hopping
            else:
                matrix[row, transmission_index] = hopping * np.exp(
                    1j * wave_number * number_of_sites
                )

        # Right lead boundary equation at j = N.
        row = number_of_sites + 1
        matrix[row, number_of_sites - 1] = hopping
        matrix[row, transmission_index] = energy * np.exp(
            1j * wave_number * number_of_sites
        ) + hopping * np.exp(1j * wave_number * (number_of_sites + 1))

        solution = np.linalg.solve(matrix, rhs)

        reflection_amplitude = solution[reflection_index]
        transmission_amplitude = solution[transmission_index]

        return ScatteringResult(
            energy=energy,
            reflection_amplitude=reflection_amplitude,
            transmission_amplitude=transmission_amplitude,
        )
