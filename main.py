from astropy.io.ascii import read
import matplotlib.pyplot as plt
import numpy as np
from opacities import calculate_opacities
from populations import calculate_populations

# Read the data
table_t5000 = read("t5000.dat", data_start=25)
# table_t8000 = read("t8000.dat" , data_start =25)

# Save the magnitudes
k = table_t5000.columns[0]
lgTauR = table_t5000.columns[1]
lgTau5 = table_t5000.columns[2]
depth = table_t5000.columns[3]
T = table_t5000.columns[4]
Pe = table_t5000.columns[5]

Ne, N_ionization, N_excitation = calculate_populations(T, Pe)

N_Hminus, N_HI, N_HII = N_ionization
N_HI_n1, N_HI_n2, N_HI_n3 = N_excitation

print("\n")
print("Para tau_ross=1 tenemos N_HI=", N_HI[40])
print("Para tau_ross=1 tenemos N_HI_n1=", N_HI_n1[40])

# Now we want to use as imput to calculate the opacities just the data orresponding to tauR = 1.0
# tauR = 1.0 -> lgTauR = 0.0 -> we want to select the row 40 of our data

T_tauR1 = T[40]
Pe_tauR1 = Pe[40]
Ne_tauR1 = Ne[40]

N_Hminus_tauR1 = N_Hminus[40]
N_HI_tauR1 = N_HI[40]
N_HII_tauR1 = N_HII[40]
N_HI_n1_tauR1 = N_HI_n1[40]
N_HI_n2_tauR1 = N_HI_n2[40]
N_HI_n3_tauR1 = N_HI_n3[40]

N_ionization_tauR1 = np.array([N_Hminus_tauR1, N_HI_tauR1, N_HII_tauR1])
N_excitation_tauR1 = np.array([N_HI_n1_tauR1, N_HI_n2_tauR1, N_HI_n3_tauR1])

k_TOTAL, k_HI, k_Hminus, k_es = calculate_opacities(
    T_tauR1, Pe_tauR1, Ne_tauR1, N_ionization_tauR1, N_excitation_tauR1
)

# print results
k_ff_HI, k_bf_HI_n1, k_bf_HI_n2, k_bf_HI_n3, k_bb_Ly_a, k_bb_Ly_b, k_bb_Bal = k_HI
k_ff_Hminus, k_bf_Hminus = k_Hminus

print("k_ff_HI", k_ff_HI.shape)
print("k_bf_HI_n3", k_bf_HI_n3.shape)
print("k_es", k_es.shape)

lamb_all = np.linspace(500, 20000, 19501)
index_1 = np.argmin(abs(lamb_all - 911.743))
index_2 = np.argmin(abs(lamb_all - 3646.973))
index_3 = np.argmin(abs(lamb_all - 8205.689))
index_4 = np.argmin(abs(lamb_all - 16444.327))

opacities = (
    k_ff_HI,
    k_bf_HI_n1,
    k_bf_HI_n2,
    k_bf_HI_n3,
    k_ff_Hminus,
    k_bf_Hminus,
)  # , k_es

# creamos las filas de la tabla de latex
for k, name in zip(
    opacities,
    [
        "$k_ff(HI)$",
        "k_bf(HI,n1)$",
        "k_bf(HI,n2)$",
        "$k_bf(HI,n3)$",
        "$k_ff(H-)$",
        "$k_bf(H-)$",
    ],
):  # , '$k_e$']):
    string = ""
    string += str(name + " & ")
    for i in [index_1, index_2, index_3, index_4]:
        string += "%.2e & " % k[i - 1]
        string += "%.2e & " % k[i + 1]
    print(string)

print("\n")
print(k_es)
