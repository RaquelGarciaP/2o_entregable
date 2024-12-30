# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 18:14:53 2024

@author: joanes


PLOTS ATMOSFERAS ESTELARES

"""

from astropy.io.ascii import read
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

from matplotlib.ticker import FuncFormatter


def format_x_labels(x, pos):
    if x == 0:
        return "0"  # Sin formato especial para el cero
    return f"{x:.1f} × 10⁷"


table_t5000 = read("t5000.dat", data_start=25)
table_t8000 = read("t8000.dat", data_start=25)

metatable = [table_t5000, table_t8000]
names = ["$T_{eff}=5000 K$", "$T_{eff}=8000 K$"]

colormap = cm.get_cmap("coolwarm", len(names))
colors = [colormap(i) for i in range(len(names))]
colors = colors[::-1]

# k = table_t5000.columns[0]
# lgTauR = table_t5000.columns[1]
# lgTau5 = table_t5000.columns[2]
# Depth = table_t5000.columns[3]
# T = table_t5000.columns[4]
# Pe = table_t5000.columns[5]
# Pg = table_t5000.columns[6]
# Prad =table_t5000.columns[7]
# Pturb = table_t5000.columns[8]

for ii in range(0, len(metatable)):
    plt.plot(
        metatable[ii]["col2"],
        metatable[ii]["col4"],
        label=f"{names[ii]}",  # Etiqueta para la leyenda
        color=colors[ii],
    )  # Asignar un color

plt.legend(fontsize=12)

# Invertir el eje x
ax = plt.gca()
ax.invert_xaxis()

# Configurar las marcas en el eje x para que aparezcan cada 2 unidades
x_ticks = np.arange(int(ax.get_xlim()[0]), int(ax.get_xlim()[1]) - 1.5, -1.5)
ax.set_xticks(x_ticks)

# Aplicar el formato personalizado al eje x
ax.xaxis.set_major_formatter(FuncFormatter(format_x_labels))


plt.ylabel(r"$\log(\tau_{\rm R})$", fontsize="12")
plt.xlabel(r"$r$", fontsize="12")

plt.gca().invert_xaxis()
plt.show()
# %%
Pe_Pg = [
    metatable[0]["col6"] / metatable[0]["col7"],
    metatable[1]["col6"] / metatable[1]["col7"],
]
Prad_Pg = [
    metatable[0]["col8"] / metatable[0]["col7"],
    metatable[1]["col8"] / metatable[1]["col7"],
]

for ii in range(0, len(metatable)):
    plt.plot(
        metatable[ii]["col2"],
        metatable[ii]["col5"],
        label=f"{names[ii]}",  # Etiqueta para la leyenda
        color=colors[ii],
    )  # Asignar un color

plt.legend(fontsize=12)
plt.ylabel(r"T[K]", fontsize="12")
plt.xlabel(r"$\log(\tau_{\rm R})$", fontsize="12")
plt.show()
# %%
for ii in range(0, len(metatable)):
    plt.plot(
        metatable[ii]["col2"],
        np.log10(metatable[ii]["col6"]),
        label=f"{names[ii]}",  # Etiqueta para la leyenda
        color=colors[ii],
    )  # Asignar un color

plt.legend(fontsize=12)
plt.ylabel(r"$\log(P_{e}[Ba])$", fontsize="12")
plt.xlabel(r"$\log(\tau_{R})$", fontsize="12")
plt.show()
"vamos a mirar si esta grafica hay que hacerlo con el log de P_{e} que si no queda todo casi al 0. también hay que reajustar el limite de ejes"
# %%
for ii in range(0, len(metatable)):
    plt.plot(
        metatable[ii]["col2"],
        Pe_Pg[ii],
        label=f"{names[ii]}",  # Etiqueta para la leyenda
        color=colors[ii],
    )  # Asignar un color

plt.legend(fontsize=12)
plt.ylabel(r"$P_{e}/P_{\rm g}$", fontsize="12")
plt.xlabel(r"$\log(\tau_{\rm R})$", fontsize="12")
plt.show()

# %%
for ii in range(0, len(metatable)):
    plt.plot(
        metatable[ii]["col2"],
        Prad_Pg[ii],
        label=f"{names[ii]}",  # Etiqueta para la leyenda
        color=colors[ii],
    )  # Asignar un color

plt.legend(fontsize=12)
plt.ylabel(r"$P_{\rm rad}/P_{\rm g}$", fontsize="12")
plt.xlabel(r"$\log(\tau_{\rm R})$", fontsize="12")
plt.show()

# %%
tau_R_t5000 = 10 ** metatable[0]["col2"]
tau_R_t8000 = 10 ** metatable[1]["col2"]

T = [
    ((3 / 4) * (5000**4) * (10 ** metatable[0]["col2"] + (2 / 3))) ** (1 / 4),
    ((3 / 4) * (8000**4) * (10 ** metatable[1]["col2"] + (2 / 3))) ** (1 / 4),
]

for ii in range(0, len(metatable)):
    plt.plot(
        metatable[ii]["col2"],
        T[ii],
        label=f"{names[ii]}",  # Etiqueta para la leyenda
        color=colors[ii],
    )  # Asignar un color


plt.legend(fontsize=12)
plt.ylabel(r"${\rm T(K)}$", fontsize="12")
plt.xlabel(r"$\log(\tau_{\rm R})$", fontsize="12")
plt.show()
