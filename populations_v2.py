import numpy as np
# import matplotlib . pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from astropy.io.ascii import read
import matplotlib.pyplot as plt

# We start reading the data at line =2
# (0 contains the labels , and 1 contains the units )

table_t5000 = read("t5000.dat" , data_start =25)
table_t8000 = read("t8000.dat" , data_start =25)

print(table_t5000.columns[1])

# Constants :
kb = 8.61 * 10**(-5)
# Fist of all I load the Data
# data = pd.read_csv('t5000.dat', delimiter='\s+', index_col=False,skiprows=6, header=0)

# I extract the data
'''k = data ['k']

lgTauR = data ['lgTauR']
lgTau5 = data ['lgTau5']
T = data ['T']
Pe = data ['Pe']'''

k = table_t5000.columns[0]

lgTauR = table_t5000.columns[1]
lgTau5 = table_t5000.columns[2]
T = table_t5000.columns[4]
Pe = table_t5000.columns[5]

Ne =[]

for x , y in zip ( Pe , T ) :
    Ne.append ( x / kb * y )
print ( Ne [40] , Ne [50])

# Soluciones de la ecuacion
C1 =[]
C2 =[]
C3 =[]
C4 =[]

# Poblaciones
NHminus =[]
NHI =[]
NHII =[]
NHeI =[]
NHeII =[]
NHeIII =[]

# Ionization energies of different elements :
XHminus =0.755
XHI =13.6
XHeI =24.6
XHeII =54.4

# Partition functions
UHminus =1
UHI =2.
UHII =1.
UHeI =1.
UHeII =2.
UHeIII =1.

# Poblaciones de Saha
def saha (N , U1 , U2 ,X , Te ) :
    return 2.07 * 10 **(-16) * N * ( U1 / U2 ) * Te **(3/2) * np.exp(X/(kb*Te))

for i in range ( len ( Ne ) ) :
    C1.append(saha( Ne[i] , UHminus, UHI, XHminus, T[i]))
    C2.append(saha( Ne[i] , UHI , UHII , XHI , T [ i ]) )
    C3.append(saha( Ne[i] , UHeI , UHeII , XHeI , T [ i ]) )
    C4.append(saha( Ne[i] , UHeII , UHeIII , XHeII , T [ i ]) )

for i in range ( len ( T ) ) :
    A = np . matrix ([[1. , - C1 [ i ] ,0. ,0. ,0. ,0.] ,[0. ,1. , - C2 [ i
    ] ,0. ,0. ,0.] ,[0. ,0. ,0. ,1. , - C3 [ i ] ,0.] ,[0. ,0. ,0. ,0. ,1. , - C4 [ i ]] ,
    [ -1. ,0. ,1. ,0. ,1. ,2.] ,[ -0.085 , -0.085 , -0.085 ,1. ,1. ,1.]])
    B = np . matrix ([[0.] ,[0.] ,[0.] ,[0.] ,[ Ne [ i ]] ,[0.]])
    Ainv = np . linalg . inv ( A )
    x = - Ainv * B
    NHminus . append ( x [0 ,0])
    NHI . append ( x [1 ,0])
    NHII . append ( x [2 ,0])
    NHeI . append ( x [3 ,0])
    NHeII . append ( x [4 ,0])
    NHeIII . append ( x [5 ,0])

N =[ NHminus , NHI , NHII , NHeI , NHeII , NHeIII ]

# Excitation energies and statistical weights for boltzmann
empty = float('NaN')
EEHI=[0 ,10.1988 ,12.0875 ,12.7485 , empty ]
gHI =[2 ,8 ,18 ,32 , empty ]
EEHeI =[0 ,19.8196 ,20.6158 ,20.9642 ,21.2180]
gHeI =[1 ,3 ,1 ,9 ,3]
EEHeII =[0 ,40.8130 , empty , empty , empty ]
gHeII =[2 ,8 , empty , empty , empty ]

EE =[ empty , EEHI , EEHeI , EEHeII ]
g =[ empty , gHI , gHeI , gHeII ]

def boltzmann (i ,j ,U , Te ) :
    n =[]
    for k in range ( len ( Te ) ) :
        n.append( N [ i ][ k ] * ( g [ i ][ j ]/ U ) * np . exp ( - EE [ i ][ j ]) / kb * Te [ k ])
        return n

nHI_1 = boltzmann (1 ,0 , UHI , T )
nHI_2 = boltzmann (1 ,1 , UHI , T )
nHI_3 = boltzmann (1 ,2 , UHI , T )
nHI_4 = boltzmann (1 ,3 , UHI , T )
nHeI_1 = boltzmann (2 ,0 , UHeI , T )
nHeI_2 = boltzmann (2 ,1 , UHeI , T )
nHeI_3 = boltzmann (2 ,2 , UHeI , T )
nHeI_4 = boltzmann (2 ,3 , UHeI , T )
nHeI_5 = boltzmann (2 ,4 , UHeI , T )
nHeII_1 = boltzmann (3 ,0 , UHeII , T )
nHeII_2 = boltzmann (3 ,1 , UHeII , T )

N = pd.DataFrame({'NHminus': NHminus , 'NHI': NHI , 'NHII': NHII , 'NHeI': NHeI, 'NHeII': NHeII, 'NHeIII': NHeIII})
n = pd.DataFrame({'nHI_1': nHI_1 , 'nHI_2': nHI_2 , 'nHI_3': nHI_3 , 'nHI_4': nHI_4 , 'nHeI_1': nHeI_1, 
                  'nHeI_2': nHeI_2 , 'nHeI_3': nHeI_3 , 'nHeI_4': nHeI_4 , 'nHeI_5': nHeI_5 , 'nHeII_1': nHeII_1,
                  'nHeII_2':nHeII_2 })
# Comprobacion
for i in range ( len ( Ne ) ) :

    print( k [ i ] , 'Conservacion de la carga' , round ( Ne [ i ]) == - round ((NHII [ i ]) - NHminus [ i ]+ NHeII [ i ]+2 * NHeIII [ i ]) )
    print('Relaciones de poblaciones He / H', round (( NHeI [ i ]+ NHeII [ i ]+ NHeIII [ i ]) /( NHminus [ i ]+ NHI [ i ]+ NHII [ i ]) ) == round (10**(10.93 -12) ) )

N = pd.read_csv('Nsaha_5000.txt', delimiter = '\s + ' , index_col = False , header =0)
n.to_csv('nboltzmann_8000.txt', header = True , index = False , sep = '\t', mode = 'a')
