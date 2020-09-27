import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### TRATAMENTO DOS DADOS DE SISMO ### 

# Pandas uses DATAFRAME objects, while Numpy uses ARRAY objects
# Is necessary convert from a object to another, using the to_numpy() method

data1 = pd.read_excel('dados_sismo.xlsx', header=None, sheet_name='comp_0').to_numpy().T
data2 = pd.read_excel('dados_sismo.xlsx', header=None, sheet_name='comp_90').to_numpy().T
data3 = pd.read_excel('dados_sismo.xlsx', header=None, sheet_name='up').to_numpy().T

def adjust(data):
    
    # Transform the seismic data into a unique column vector
    
    lin, col = np.shape(data)
    fited = np.zeros((lin*col,1))
    
    for i in range(0,col):
        aux1 = 8*i
        aux2 = aux1+lin
        fited[aux1:aux2,0] = np.copy(data[0:lin,i])

    return fited

comp_0 = adjust(data1)
comp_90 = adjust(data2)
up = adjust(data3)

seismic = np.zeros((np.shape(comp_0)[0],3))
seismic[:,0] = np.copy(comp_0[:,0])
seismic[:,1] = np.copy(comp_90[:,0])
seismic[:,2] = np.copy(up[:,0])

real_seismic = seismic.T

df     = pd.DataFrame(real_seismic)
writer = pd.ExcelWriter('sismo_real.xlsx')
df.to_excel(writer,'Sheet1', index=False) 
writer.save()

t = np.linspace(0,39.98,num=2000)

### PLOTTING GRAPHICS ###

plt.figure(1, figsize=(12,3)) 
plt.plot(t,0.01*comp_0,'b-')
plt.title("Componente 0°")
plt.xlabel('tempo (s)')
plt.ylabel('Aceleração (m/s²)')
plt.xlim(0,39.98);plt.ylim(-0.3,0.3)
plt.grid(True)

plt.figure(2, figsize=(12,3)) 
plt.plot(t,0.01*comp_90,'r')
plt.title("Componente 90°")
plt.xlabel('tempo (s)')
plt.ylabel('Aceleração (m/s²)')
plt.xlim(0,39.98);plt.ylim(-0.7,0.6)
plt.grid(True)

plt.figure(3, figsize=(12,3)) 
plt.plot(t,0.01*up,'g')
plt.title("Componente UP")
plt.xlabel('tempo (s)')
plt.ylabel('Aceleração (m/s²)')
plt.xlim(0,39.98);plt.ylim(-0.3,0.3)
plt.grid(True)

plt.show()
