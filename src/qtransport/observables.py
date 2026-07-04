def reflection_probability(result):
    return abs(result["r"]) ** 2


def transmission_probability(result):
    return abs(result["t"]) ** 2


def conservation_error(result):
    R = reflection_probability(result)
    T = transmission_probability(result)
    return abs(R + T - 1)
