import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import pandas as pd
import Dynas as dyn 

###    ENTRADA DE DADOS     ###

t = np.linspace(0,39.98,num=2000)

u_a = pd.read_excel('deslocamentos_artificial.xlsx').to_numpy()
u_r = pd.read_excel('deslocamentos_real.xlsx').to_numpy()

### MODELANDO O DRIFT PARA SISMO ARTIFICIAL ###

drift_u_a = np.zeros((8,2000))
drift_u_a[0,:] = u_a[138,:]

for i in range (0,2000):
    for k in range (0,7):
        drift_u_a[k+1,i] = u_a[(138+144*(k+1)),i] - u_a[(138+144*k),i]  

for i in range (0,8):
    plt.figure(i,figsize=(12,4))   
    plt.plot(t,100*drift_u_a[i,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Deslocamento relativo entre pavimentos (cm)')
    plt.xlim(0,max(t)); plt.ylim(-0.20,0,20)
    plt.title(f'Deslocamento relativo entre o {i}º e o {i+1}º pavimentos na direção 0 DEG')
    plt.grid(True)
    plt.show()
