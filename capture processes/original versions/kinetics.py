import math

# References
T_REF = 313.15  # more realistic reference for carbon capture co2 column
R = 8.314462618

# Pre-exponential factor
A = 4.0e6
E_A = 41000.0  # J/mol


def arr_constant(t_user):
    """Calculate the temperature-dependent reaction rate constant using the

    centered Arrhenius equation relative to T_REF:

        k = A * exp(-(Ea / R) * (1/T - 1/T_REF))
    """

    if t_user <= 0:
        raise ValueError("Temperature must be greater than zero.")

    # FIXED: Incorporated T_REF into the exponential term
    constant = A * math.exp(-(E_A / R) * ((1.0 / t_user) - (1.0 / T_REF)))
    return constant


# CO2 + MEA REACTION RATE


def reaction_rate(t_user, co2_conc, mea_conc):
    """Calculate the reaction rate between dissolved CO2 and free MEA.

        r = k [CO2][MEA]
    """
    if co2_conc < 0:
        raise ValueError("CO2 concentration cannot be negative.")

    if mea_conc < 0:
        raise ValueError("MEA concentration cannot be negative.")

    k = arr_constant(t_user)

    rate = k * co2_conc * mea_conc
    return rate
