import numpy as np
import numpy.linalg as la
import scipy.linalg as sc
import matplotlib.pyplot as plt
import pandas as pd
from pandas import ExcelWriter
from numba import jit
from mpl_toolkits.mplot3d import Axes3D 
import Dynas as dyn 

# Lista de graus de liberdade restringidos

Restrictions = np.arange(144)

# Listando as frequencias naturais

K, M = dyn.Matrix3D('dados.xlsx',216)
Kr, Mr = dyn.Restr(K, M, Restrictions)
fk, wk, Phi = dyn.modal_analysis(Kr, Mr, 1152)

# C = dyn.damping_matrix(wk[0], wk[1], 0.01, 0.01, M, K)

### PLOTANDO O GRÁFICO ###

"""plt.subplot(4,1,1)
plt.plot(t,df[0,:],'orange') # plotando o deslocamento do primeiro grau de liberdade (73º GL)
plt.title('Diferenças finitas centrais')
plt.grid(True)"""