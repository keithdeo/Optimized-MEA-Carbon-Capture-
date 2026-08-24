import math

from properties import R


# KINETIC PARAMETERS

# Reference temperature
T_REF = 313.15       # K


# Pre-exponential factor
A = 4.0e6


# Activation energy
E_A = 41000.0        # J/mol


# ARRHENIUS RATE CONSTANT

def rate_constant(
    temperature_K
):
    """
    Calculate the temperature-dependent reaction
    rate constant using the Arrhenius equation:

        k = A exp(-Ea / RT)
    """

    if temperature_K <= 0:

        raise ValueError(
            "Temperature must be greater than zero."
        )

    return (
        A
        * math.exp(
            -E_A
            /
            (R * temperature_K)
        )
    )


# CO2 + MEA REACTION RATE

def reaction_rate(
    co2_concentration,
    mea_concentration,
    temperature_K
):
    """
    Calculate the reaction rate between dissolved
    CO2 and free MEA.

        r = k [CO2][MEA]
    """

    if co2_concentration < 0:

        raise ValueError(
            "CO2 concentration cannot be negative."
        )

    if mea_concentration < 0:

        raise ValueError(
            "MEA concentration cannot be negative."
        )

    k = rate_constant(
        temperature_K
    )

    return (
        k
        * co2_concentration
        * mea_concentration
    )
