from astropy.io.ascii import read
import matplotlib.pyplot as plt
import numpy as np
from opacities import calculate_opacities
from populations import calculate_populations

# from prueba import calculate_populations

# *********** Read the data ***********
table_data = read("t5000.dat", data_start=25)
# table_data = read("t8000.dat" , data_start =25)

# Save the magnitudes
k = table_data.columns[0]
lgTauR = table_data.columns[1]
lgTau5 = table_data.columns[2]
depth = table_data.columns[3]
T = table_data.columns[4]
Pe = table_data.columns[5]

# *********** Calculate populations ***********

Ne, N_ionization, N_excitation = calculate_populations(T, Pe)

N_Hminus, N_HI, N_HII = N_ionization
N_HI_n1, N_HI_n2, N_HI_n3 = N_excitation


# Now we want to use as imput to calculate the opacities just the data orresponding to tauR = 1.0
# tauR = 1.0 -> lgTauR = 0.0 -> we want to select the row 40 of our data

index_tauR1 = 40

T_tauR1 = T[index_tauR1]
Pe_tauR1 = Pe[index_tauR1]
Ne_tauR1 = Ne[index_tauR1]

N_Hminus_tauR1 = N_Hminus[index_tauR1]
N_HI_tauR1 = N_HI[index_tauR1]
N_HII_tauR1 = N_HII[index_tauR1]
N_HI_n1_tauR1 = N_HI_n1[index_tauR1]
N_HI_n2_tauR1 = N_HI_n2[index_tauR1]
N_HI_n3_tauR1 = N_HI_n3[index_tauR1]

N_ionization_tauR1 = np.array([N_Hminus_tauR1, N_HI_tauR1, N_HII_tauR1])
N_excitation_tauR1 = np.array([N_HI_n1_tauR1, N_HI_n2_tauR1, N_HI_n3_tauR1])

# *********** Calculate opacities ***********

k_TOTAL, k_HI, k_Hminus, k_es = calculate_opacities(
    T_tauR1, Pe_tauR1, Ne_tauR1, N_ionization_tauR1, N_excitation_tauR1
)


# *********** Print all the results ***********

# Populations
print("*********** Results populations ***********")

index_tauR05 = 37
index_tauR5 = 47

N_populations = N_Hminus, N_HI, N_HII, Ne, N_HI_n1, N_HI_n2, N_HI_n3

# we create the rows of a latex table with the results
print("\n Table with populations for TauR = 0.5 and TauR = 5 \n")

print("$\ tau$ & H- & HI & HII & Ne & HI, n=1 & HI, n=2 & HI, n=3")
for i, tau in zip([index_tauR05, index_tauR5], ["$0.5$", "$5$"]):
    string = ""
    string += str(tau + " & ")
    for population in N_populations:
        string += "%.2e & " % population[i]
    print(string)


# Opacities
print("\n")
print("*********** Results opacities ***********")
k_ff_HI, k_bf_HI_n1, k_bf_HI_n2, k_bf_HI_n3, k_bb_Ly_a, k_bb_Ly_b, k_bb_Bal = k_HI
k_ff_Hminus, k_bf_Hminus = k_Hminus

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
)

# we create the rows of a latex table with the results
print("\n Table with opacities k_ff and k_bf of HI and H- \n")
for k, name in zip(
    opacities,
    [
        "$k_{ff}(HI)$",
        "k_{bf}(HI,n1)$",
        "k_{bf}(HI,n2)$",
        "$k_{bf}(HI,n3)$",
        "$k_{ff}(H-)$",
        "$k_{bf}(H-)$",
    ],
):  # , '$k_e$']):
    string = ""
    string += str(name + " & ")
    for i in [index_1, index_2, index_3, index_4]:
        string += "%.2e & " % k[i - 1]
        string += "%.2e & " % k[i + 1]
    print(string)

print("$k_e$ & %.2e" % k_es)

print("\nLambdas used:")
print("Lamda_1 - Delta = ", lamb_all[index_1 - 1])
print("Lamda_1 + Delta = ", lamb_all[index_1 + 1])
print("Lamda_2 - Delta = ", lamb_all[index_2 - 1])
print("Lamda_2 + Delta = ", lamb_all[index_2 + 1])
print("Lamda_3 - Delta = ", lamb_all[index_3 - 1])
print("Lamda_3 + Delta = ", lamb_all[index_3 + 1])
print("Lamda_4 - Delta = ", lamb_all[index_4 - 1])
print("Lamda_4 + Delta = ", lamb_all[index_4 + 1])

# now we print the table of opacities for Lyman and Balmer
opacities_lyman_balmer = [k_bb_Ly_a, k_bb_Ly_b, k_bb_Bal]
print("\n Table with opacities for Lyman alpha, Lyman beta and Balmer alpha \n")

print(" & $L_{\ alpha}$ & $L_{\ beta}$ & $H_{\ alpha}$")
string = " & %.2e " % k_bb_Ly_a
string += "& %.2e " % k_bb_Ly_b
string += "& %.2e " % k_bb_Bal
print(string)
