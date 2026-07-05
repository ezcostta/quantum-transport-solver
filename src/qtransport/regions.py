import numpy as np


class OneDimensionalRegion:
    def __init__(self, onsite, hopping: float = 1.0):
        self.onsite = np.array(onsite, dtype=float)
        self.hopping = hopping

    @property
    def size(self) -> int:
        return len(self.onsite)
