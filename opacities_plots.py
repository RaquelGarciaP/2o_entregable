from astropy.io.ascii import read
import matplotlib.pyplot as plt
import numpy as np
from opacities import calculate_opacities
from populations import calculate_populations

# *********** Read the data ***********
table_t5000 = read("t5000.dat", data_start=25)
table_t8000 = read("t8000.dat", data_start=25)

# Save the magnitudes for table_t5000
T_5000 = table_t5000.columns[4]
Pe_5000 = table_t5000.columns[5]

# Save the magnitudes for table_t8000
T_8000 = table_t8000.columns[4]
Pe_8000 = table_t8000.columns[5]

# *********** Calculate populations ***********

# For table_t5000
Ne_5000, N_ionization_5000, N_excitation_5000 = calculate_populations(T_5000, Pe_5000)

# For table_t8000
Ne_8000, N_ionization_8000, N_excitation_8000 = calculate_populations(T_8000, Pe_8000)

# *********** Calculate opacities for tauR = 1.0 ***********

# Now we want to use as imput to calculate the opacities just the data orresponding to tauR = 1.0
# tauR = 1.0 -> lgTauR = 0.0 -> we want to select the row 40 of our data
index_tauR1 = 40

# Calculate opacities for table_t5000
T_tauR1_5000 = T_5000[index_tauR1]
Pe_tauR1_5000 = Pe_5000[index_tauR1]
Ne_tauR1_5000 = Ne_5000[index_tauR1]

N_ionization_5000 = np.array(N_ionization_5000)
N_excitation_5000 = np.array(N_excitation_5000)

N_ionization_tauR1_5000 = N_ionization_5000[:, index_tauR1]
N_excitation_tauR1_5000 = N_excitation_5000[:, index_tauR1]

k_TOTAL_5000, k_HI_5000, k_Hminus_5000, k_es_5000 = calculate_opacities(
    T_tauR1_5000,
    Pe_tauR1_5000,
    Ne_tauR1_5000,
    N_ionization_tauR1_5000,
    N_excitation_tauR1_5000,
)

(
    k_ff_HI_5000,
    k_bf_HI_n1_5000,
    k_bf_HI_n2_5000,
    k_bf_HI_n3_5000,
    k_bb_Ly_a_5000,
    k_bb_Ly_b_5000,
    k_bb_Bal_5000,
) = k_HI_5000
k_ff_Hminus_5000, k_bf_Hminus_5000 = k_Hminus_5000

# Calculate opacities for table_t8000
T_tauR1_8000 = T_8000[index_tauR1]
Pe_tauR1_8000 = Pe_8000[index_tauR1]
Ne_tauR1_8000 = Ne_8000[index_tauR1]

N_ionization_8000 = np.array(N_ionization_8000)
N_excitation_8000 = np.array(N_excitation_8000)

N_ionization_tauR1_8000 = N_ionization_8000[:, index_tauR1]
N_excitation_tauR1_8000 = N_excitation_8000[:, index_tauR1]

k_TOTAL_8000, k_HI_8000, k_Hminus_8000, k_es_8000 = calculate_opacities(
    T_tauR1_8000,
    Pe_tauR1_8000,
    Ne_tauR1_8000,
    N_ionization_tauR1_8000,
    N_excitation_tauR1_8000,
)

(
    k_ff_HI_8000,
    k_bf_HI_n1_8000,
    k_bf_HI_n2_8000,
    k_bf_HI_n3_8000,
    k_bb_Ly_a_8000,
    k_bb_Ly_b_8000,
    k_bb_Bal_8000,
) = k_HI_8000
k_ff_Hminus_8000, k_bf_Hminus_8000 = k_Hminus_8000

# *********** Plots ***********

lamb_all = np.linspace(500, 20000, 19501)

# Define font size of the plots
plt.rcParams.update(
    {
        "font.size": 18,
        "axes.titlesize": 16,
        "axes.labelsize": 18,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
    }
)

# Plot opacities for T = 5000 K
plt.figure(figsize=(12, 8))
plt.plot(
    lamb_all,
    k_TOTAL_5000,
    label=r"$k_{\rm TOTAL}$",
    linestyle="-",
    color="black",
    linewidth=2,  # make k_TOTAL bigger
)
plt.plot(
    lamb_all,
    k_ff_HI_5000,
    label=r"$k_{ff}$ $({\rm HI})$",
    linestyle=":",
    color="darkturquoise",
)
plt.plot(
    lamb_all,
    k_bf_HI_n1_5000,
    label=r"$k_{bf}$ $({\rm HI}, n=1)$",
    linestyle=":",
    color="yellowgreen",
)
plt.plot(
    lamb_all,
    k_bf_HI_n2_5000,
    label=r"$k_{bf}$ $({\rm HI}, n=2)$",
    linestyle="-.",
    color="darkorange",
)
plt.plot(
    lamb_all,
    k_bf_HI_n3_5000,
    label=r"$k_{bf}$ $({\rm HI}, n=3)$",
    linestyle=":",
    color="darkred",
)
plt.plot(
    lamb_all,
    k_ff_Hminus_5000,
    label=r"$k_{ff}$ $({\rm H}^-)$",
    linestyle="--",
    color="purple",
)
plt.plot(
    lamb_all,
    k_bf_Hminus_5000,
    label=r"$k_{bf}$ $({\rm H}^-)$",
    linestyle="--",
    color="green",
)
plt.plot(
    lamb_all,
    np.full_like(lamb_all, k_es_5000),
    label=r"$k_{e}$",
    color="navy",
    linestyle="-.",
)

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\lambda$ $(\AA)$")
plt.ylabel(r"$k$ ${\rm (cm^{-1})}$")
plt.legend()
# plt.grid(which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()


# Plot opacities for T = 8000 K
plt.figure(figsize=(12, 8))
plt.plot(
    lamb_all,
    k_TOTAL_8000,
    label=r"$k_{\rm TOTAL}$",
    linestyle="-",
    color="black",
    linewidth=2,  # make k_TOTAL bigger
)
plt.plot(
    lamb_all,
    k_ff_HI_8000,
    label=r"$k_{ff}$ $({\rm HI})$",
    linestyle=":",
    color="darkturquoise",
)
plt.plot(
    lamb_all,
    k_bf_HI_n1_8000,
    label=r"$k_{bf}$ $({\rm HI}, n=1)$",
    linestyle=":",
    color="yellowgreen",
)
plt.plot(
    lamb_all,
    k_bf_HI_n2_8000,
    label=r"$k_{bf}$ $({\rm HI}, n=2)$",
    linestyle="-.",
    color="darkorange",
)
plt.plot(
    lamb_all,
    k_bf_HI_n3_8000,
    label=r"$k_{bf}$ $({\rm HI}, n=3)$",
    linestyle=":",
    color="darkred",
)
plt.plot(
    lamb_all,
    k_ff_Hminus_8000,
    label=r"$k_{ff}$ $({\rm H}^-)$",
    linestyle="--",
    color="purple",
)
plt.plot(
    lamb_all,
    k_bf_Hminus_8000,
    label=r"$k_{bf}$ $({\rm H}^-)$",
    linestyle="--",
    color="green",
)
plt.plot(
    lamb_all,
    np.full_like(lamb_all, k_es_8000),
    label=r"$k_{e}$",
    color="navy",
    linestyle="-.",
)

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\lambda$ $(\AA)$")
plt.ylabel(r"$k$ ${\rm (cm^{-1})}$")
plt.legend()
# plt.grid(which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()


# *********** Plot total opacity vs wavelength ***********

plt.figure(figsize=(10, 6))
plt.plot(lamb_all, k_TOTAL_5000, label="T = 5000 K", color="red")
plt.plot(lamb_all, k_TOTAL_8000, label="T = 8000 K", color="blue")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$\lambda$ $(\AA)$")
plt.ylabel(r"$k$ ${\rm (cm^{-1})}$")
plt.legend()
# plt.grid(which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()
# plt.close()
