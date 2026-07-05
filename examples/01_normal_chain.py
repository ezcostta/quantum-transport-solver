from qtransport.solver import OneDimensionalScatteringSolver
from qtransport.observables import (
    reflection_probability,
    transmission_probability,
    conservation_error,
)


def main():
    onsite_region = [0.5, 0.5, 0.5, 0.5]
    solver = OneDimensionalScatteringSolver(
        onsite_region=onsite_region,
        hopping=1.0,
    )

    energy = 0.3
    result = solver.solve(energy)

    R = reflection_probability(result)
    T = transmission_probability(result)

    print(f"Energy       = {result['energy']:.6f}")
    print(f"k            = {result['k']:.6f}")
    print(f"r            = {result['r']:.6f}")
    print(f"t            = {result['t']:.6f}")
    print(f"R            = {R:.12f}")
    print(f"T            = {T:.12f}")
    print(f"R + T        = {R + T:.12f}")
    print(f"error        = {conservation_error(result):.2e}")


if __name__ == "__main__":
    main()
