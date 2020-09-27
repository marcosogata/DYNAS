import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import pandas as pd
import Dynas as dyn 

# Lista de graus de liberdade restringidos

Restrictions = np.arange(144)

# Montando as Matrizes de Rigidez e de Massa

K, M = dyn.Matrix3D('dados.xlsx',216)
Kr, Mr = dyn.Restr(K, M, Restrictions)

fk, wk, Phi = dyn.modal_analysis(Kr, Mr, 4)

#################################
###      CARGAS DE SISMO      ###
#################################

# Os arquivos de entrada de dados estão em cm/s/s. É necessário converter pra SI /100

# Obtendo a carga de sismo real

sismo_real = pd.read_excel('sismo_real.xlsx').to_numpy()/100
t = np.linspace(0,39.98,num=2000)
F = dyn.Seismic3D('forca_sismo_real',Mr,sismo_real,t)

# Obtendo a carga de sismo artificial

comp_0 = pd.read_excel('sismo_artificial_comp0.xlsx').to_numpy().T
comp_90 = pd.read_excel('sismo_artificial_comp90.xlsx').to_numpy().T
up = pd.read_excel('sismo_artificial_up.xlsx').to_numpy().T

col = np.shape(comp_0)[1]
sismo_artificial = np.zeros((3,col))

sismo_artificial[0,:] = np.copy(comp_0)
sismo_artificial[1,:] = np.copy(comp_90)
sismo_artificial[2,:] = np.copy(up)

sismo_artificial /= 100

F = dyn.Seismic3D('forca_sismo_artificial',Mr,sismo_artificial,t)