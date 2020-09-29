import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sc

t = np.linspace(0,39.98,num=2000)

a_sismo_real = pd.read_excel('sismo_real.xlsx').to_numpy()

a_sismo_artificial_0 = pd.read_excel('sismo_artificial_comp0.xlsx').to_numpy().T
a_sismo_artificial_90 = pd.read_excel('sismo_artificial_comp90.xlsx').to_numpy().T
a_sismo_artificial_up = pd.read_excel('sismo_artificial_up.xlsx').to_numpy().T

plt.figure(1,figsize=(12,4))
f,Saz = sc.periodogram(a_sismo_real[0,:],(len(a_sismo_real[0,:])/np.max(t)))
plt.plot(f,Saz,'r')
plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)')
plt.xlim(0,20); plt.title(' Espectro de aceleração do sismo real na direção 0DEG')
plt.grid(True)
plt.show()

plt.figure(2,figsize=(12,4))
f,Saz = sc.periodogram(a_sismo_real[1,:],(len(a_sismo_real[1,:])/np.max(t)))
plt.plot(f,Saz,'r')
plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)')
plt.xlim(0,20); plt.title(' Espectro de aceleração do sismo real na direção 90DEG')
plt.grid(True)
plt.show()

plt.figure(3,figsize=(12,4))
f,Saz = sc.periodogram(a_sismo_real[2,:],(len(a_sismo_real[2,:])/np.max(t)))
plt.plot(f,Saz,'r')
plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)')
plt.xlim(0,20); plt.title(' Espectro de aceleração do sismo real na direção UP')
plt.grid(True)
plt.show()

plt.figure(4,figsize=(12,4))
f,Saz = sc.periodogram(a_sismo_artificial_0[0,:],(len(a_sismo_artificial_0[0,:])/np.max(t)))
plt.plot(f,Saz,'r')
plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)')
plt.xlim(0,20); plt.title(' Espectro de aceleração do sismo artificial na direção 0DEG')
plt.grid(True)
plt.show()

plt.figure(5,figsize=(12,4))
f,Saz = sc.periodogram(a_sismo_artificial_90[0,:],(len(a_sismo_artificial_90[0,:])/np.max(t)))
plt.plot(f,Saz,'r')
plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)')
plt.xlim(0,20); plt.title(' Espectro de aceleração do sismo artificial na direção 90DEG')
plt.grid(True)
plt.show()

plt.figure(6,figsize=(12,4))
f,Saz = sc.periodogram(a_sismo_artificial_up[0,:],(len(a_sismo_artificial_up[0,:])/np.max(t)))
plt.plot(f,Saz,'r')
plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)')
plt.xlim(0,20); plt.title(' Espectro de aceleração do sismo artificial na direção UP')
plt.grid(True)
plt.show()