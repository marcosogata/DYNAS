import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import pandas as pd
import Dynas as dyn 

###    ENTRADA DE DADOS     ###

t = np.linspace(0,39.98,num=2000)

Mr = pd.read_excel('matriz de massa.xlsx').to_numpy()
Kr = pd.read_excel('matriz de rigidez.xlsx').to_numpy()

u0 = np.zeros((np.shape(Mr)[0],1))
v0 = np.zeros((np.shape(Mr)[0],1))

###      ANÁLISE MODAL      ###

fk, wk, Phi = dyn.modal_analysis(Kr, Mr, 4)

### MATRIZ DE AMORTECIMENTO ###

Cr = dyn.damping_matrix(wk[1], wk[3], 0.01, 0.01, Mr, Kr)

##################################################
###   RESOLUÇÃO DO PROBLEMA PARA SISMO REAL    ###
##################################################

F_real = pd.read_excel('forca_sismo_real.xlsx').to_numpy()

u_real, v_real, a_real = dyn.Newmark('real',Kr,Mr,Cr,F_real,u0,v0,t)

###################################################
### RESOLUÇÃO DO PROBLEMA PARA SISMO ARTIFICIAL ###
###################################################

F_artificial = pd.read_excel('forca_sismo_artificial.xlsx').to_numpy()

u_artificial, v_artificial, a_artificial = dyn.Newmark('artificial',Kr,Mr,Cr,F_artificial,u0,v0,t)
