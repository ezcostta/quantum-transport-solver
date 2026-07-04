import numpy as np


class OneDimensionalLead:
    """
    Semi-infinite 1D tight-binding lead with dispersion

        E = eps - 2 t cos(k)

    Hamiltonian convention:
        H = eps c†c - t(c†_{j+1}c_j + h.c.)
    """

    def __init__(self, hopping: float = 1.0, onsite: float = 0.0):
        self.hopping = hopping
        self.onsite = onsite

    def wave_number(self, energy: float) -> float:
        argument = (self.onsite - energy) / (2 * self.hopping)

        if abs(argument) > 1:
            raise ValueError("Energy is outside the propagating band.")

        return np.arccos(argument)

    def velocity(self, energy: float) -> float:
        k = self.wave_number(energy)
        return 2 * self.hopping * np.sin(k)
