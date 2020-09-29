import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

t = np.linspace(0,39.98,num=2000)
u = pd.read_excel('deslocamentos_artificial.xlsx').to_numpy()

def plot_graph(u,t,no):
    """u  = matrix of displacements, accelerations or velocitys
       no = number of node"""

    plt.figure(1,figsize=(12,4))   
    plt.plot(t,100*u[no*6-144-6 ,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Deslocamento (cm)')
    plt.xlim(0,max(t))
    plt.ylim(-max(100*u[no*6-144-6,:])*1.1,max(100*u[no*6-144-6,:])*1.1)
    plt.title(f'Deslocamento do nó {no} na direção 0 DEG')
    plt.grid(True)
    plt.show()

    plt.figure(2,figsize=(12,4))   
    plt.plot(t,100*u[no*6-144-5 ,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Deslocamento (cm)')
    plt.xlim(0,max(t))
    plt.ylim(-max(100*u[no*6-144-5,:])*1.1,max(100*u[no*6-144-5,:])*1.1)
    plt.title(f'Deslocamento do nó {no} na direção 90 DEG')
    plt.grid(True)
    plt.show()

    plt.figure(3,figsize=(12,4))   
    plt.plot(t,100*u[no*6-144-4 ,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Deslocamento (cm)')
    plt.xlim(0,max(t))
    plt.ylim(-max(100*u[no*6-144-4,:])*1.1,max(100*u[no*6-144-4,:])*1.1)
    plt.title(f'Deslocamento do nó {no} na direção UP')
    plt.grid(True)
    plt.show()

plot_graph(u,t,212)
plot_graph(u,t,201)
plot_graph(u,t,195)
plot_graph(u,t,213)
plot_graph(u,t,202)