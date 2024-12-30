import matplotlib.pyplot as plt
import numpy as np

# Constants
kb_eV = 8.6173324e-5  # Boltzmann constant (eV*K^-1)
kb = 1.380649e-16  # erg/K = g*cm^2/s*K
h = 6.6261e-27  # erg*s = cm^2*g/s
c = 2.998e10  # cm/s
e = 4.803e-10  # electron charge (esu = cm^3/2 * g^1/2 * s^-1)
m = 9.1094e-28  # electron mass (g)

sigma_e = 6.25e-25  # electron cross-section (cm^2)
R = 1.0968e5  # Rydberg constant (cm^-1)

# complete wavelength array (from 500 to 20000 A)
lamb_A = np.linspace(500, 20000, 19501)  # A (19501 so we have a step of 1 A)
lamb = lamb_A * 1e-8  # cm
# frequency array
nu = c / lamb

# a_i constants of the H- bound-free absorption cross-section
aa = [
    1.99654,
    -1.18267e-5,
    2.64243e-6,
    -4.40524e-10,
    3.23992e-14,
    -1.39568e-18,
    2.78701e-23,
]


# function to calculate the opacities
def calculate_opacities(T, Pe, Ne, N_ionization, N_excitation):

    # Save the ionization and excitation populations
    N_Hminus, N_HI, N_HII = N_ionization
    N_HI_n1, N_HI_n2, N_HI_n3 = N_excitation

    # *********** H (HI) opacities ***********

    # free-free
    g_ff = 1 + (0.3456 / (lamb * R) ** 3) * ((lamb * kb * T) / (h * c) + 1 / 2)
    sigma_ff_HI = 3.7e8 * g_ff / (np.sqrt(T) * nu**3)

    k_ff_HI = sigma_ff_HI * Ne * N_HI * (1 - np.exp(-h * nu / (kb * T)))

    # bound-free
    nn = [1, 2, 3]
    k_bf_HI_nn = []  # array where we will save the opacity for each excitation level

    # loop to calculate the opacities for HI n=1,2,3
    for n, N_HI_nn in zip(nn, N_excitation):

        # Condition for lamb values
        lamb_condition = lamb <= (n**2 / R)  # Boolean mask for lamb condition

        # Initialize k_bf_HI to zeros
        k_bf_HI = np.zeros_like(lamb)

        # Calculate for lamb values that satisfy the condition
        g_bf = 1 + (0.3456 / (lamb[lamb_condition] * R) ** 3) * (
            (lamb[lamb_condition] * R) / n**2 - 1 / 2
        )
        sigma_bf_HI = 2.815e29 * g_bf / (n**5 * nu[lamb_condition] ** 3)

        # Update k_bf_HI only where condition is met
        k_bf_HI[lamb_condition] = (
            sigma_bf_HI * N_HI_nn * (1 - np.exp(-h * nu[lamb_condition] / (kb * T)))
        )

        # Append the calculated k_bf_HI for this n
        k_bf_HI_nn.append(k_bf_HI)

    # bound-bound (line absorption -> Balmer and Lyman series)

    def line_absorption_coef(u, l, n_u, n_l, delta):
        g_bb = (
            np.pi
            * np.sqrt(3)
            * abs((((u - l) / (u + l)) ** (2 * u + 2 * l) * u * l * delta) / (u - l))
        )
        f = (
            2**5
            * g_bb
            / (3 ** (3 / 2) * np.pi * l**5 * u**3)
            * (1 / l**2 - 1 / u**2) ** (-3)
        )
        sigma_bb = np.pi * e**2 * f / (m * c)

        gu = 2 * u**2
        gl = 2 * l**2

        k_bb = sigma_bb * (n_l - gl / gu * n_u)

        return k_bb

    # Lyman alpha
    u = 2
    l = 1
    delta_Ly_a = 8 * u * (u + 1) / (u - 1) ** 2

    k_bb_Ly_a = line_absorption_coef(u, l, N_HI_n2, N_HI_n1, delta_Ly_a)

    # Lyman beta
    u = 3
    l = 1
    delta_Ly_b = 8 * u * (u + 1) / (u - 1) ** 2

    k_bb_Ly_b = line_absorption_coef(u, l, N_HI_n3, N_HI_n1, delta_Ly_b)

    # Balmer (H alpha)
    u = 3
    l = 2
    delta_Bal = 16 * u * (u + 2) * (3 * u**2 - 4) * (5 * u**2 - 4) / (u - 2) ** 6

    k_bb_Bal = line_absorption_coef(u, l, N_HI_n3, N_HI_n2, delta_Bal)

    # *********** H- opacities ***********
    theta = 5040 / T

    # free-free
    f0 = (
        -2.2763
        - 1.6850 * np.log10(lamb_A)
        + 0.76661 * ((np.log10(lamb_A)) ** 2)
        - 0.053346 * ((np.log10(lamb_A)) ** 3)
    )
    f1 = (
        15.2827
        - 9.2846 * np.log10(lamb_A)
        + 1.99381 * ((np.log10(lamb_A)) ** 2)
        - 0.142631 * ((np.log10(lamb_A)) ** 3)
    )
    f2 = (
        -197.789
        + 190.266 * np.log10(lamb_A)
        - 67.9775 * ((np.log10(lamb_A)) ** 2)
        + 10.6913 * ((np.log10(lamb_A)) ** 3)
        - 0.625151 * ((np.log10(lamb_A)) ** 4)
    )

    sigma_ff_Hminus = 1e-26 * 10 ** (
        f0 + f1 * np.log10(theta) * f2 * (np.log10(theta)) ** 2
    )

    k_ff_Hminus = sigma_ff_Hminus * Pe * N_HI

    # bound-free

    # Condition for lamb values
    e_ionization_Hminus = 1.209643358e-12  # ergs; ionization energy of H- (=0.755 eV)
    lamb_ionization_Hminus = (
        h * c / e_ionization_Hminus
    )  # different condition from HI because it is
    # not an hydrogenoid atom (has 2 e-)
    lamb_condition = lamb <= lamb_ionization_Hminus  # Boolean mask for lamb condition

    # Initialize k_bf_Hminus to zeros
    k_bf_Hminus = np.zeros_like(lamb)

    sigma_bf_Hminus = (
        aa[0]
        + aa[1] * lamb_A[lamb_condition]
        + aa[2] * lamb_A[lamb_condition] ** 2
        + aa[3] * lamb_A[lamb_condition] ** 3
        + aa[4] * lamb_A[lamb_condition] ** 4
        + aa[5] * lamb_A[lamb_condition] ** 5
        + aa[6] * lamb_A[lamb_condition] ** 6
    ) * 1e-18

    # Update k_bf_Hminus only where condition is met
    k_bf_Hminus[lamb_condition] = (
        4.158e-10 * sigma_bf_Hminus * Pe * theta ** (5 / 2) * 10 ** (0.754 * theta)
    )

    # *********** Electron scattering opacity ***********
    sigma_es = 6.25e-25
    k_es = sigma_es * Ne

    # *********** TOTAL absorption (sum of all the absorptions) ***********
    k_bf_HI_n1, k_bf_HI_n2, k_bf_HI_n3 = k_bf_HI_nn

    # save the opacities in arrays for each element (H-, HI and e-)
    k_HI = [k_ff_HI, k_bf_HI_n1, k_bf_HI_n2, k_bf_HI_n3, k_bb_Ly_a, k_bb_Ly_b, k_bb_Bal]
    k_Hminus = [k_ff_Hminus, k_bf_Hminus]

    # calculate the total opacity (sum of all opacities)
    k_TOTAL = (
        k_ff_HI
        + k_bf_HI_n1
        + k_bf_HI_n2
        + k_bf_HI_n3
        + k_ff_Hminus
        + k_bf_Hminus
        + k_es
    )

    return k_TOTAL, k_HI, k_Hminus, k_es
