from astropy.io.ascii import read
import matplotlib.pyplot as plt
import numpy as np
from populations import calculate_populations

# Read the data
table_t5000 = read("t5000.dat" , data_start =25)
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

print('\n')
print('Para tau_ross=1 tenemos N_HI=', N_HI[40])
print('Para tau_ross=1 tenemos N_HI_n1=', N_HI_n1[40])