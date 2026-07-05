from dataclasses import dataclass


@dataclass(frozen=True)
class ScatteringResult:
    energy: float
    reflection_amplitude: complex
    transmission_amplitude: complex

    @property
    def reflection(self) -> float:
        return abs(self.reflection_amplitude) ** 2

    @property
    def transmission(self) -> float:
        return abs(self.transmission_amplitude) ** 2

    @property
    def conservation_error(self) -> float:
        return abs(self.reflection + self.transmission - 1.0)
