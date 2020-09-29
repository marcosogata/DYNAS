import numpy as np
import numpy.linalg as la  
import matplotlib.pyplot as plt
import pandas as pd

F0 = np.zeros((4,3))
F0[0,2] = 10
F0[3,2] = 6
Fa = np.c_[F0, np.zeros((4,6-2))] #c_ acrescenta coluna
Fa = np.r_[Fa, [np.ones(np.shape(Fa)[1])]] #r_ acrescenta linha
print(np.shape(F0)[0])
print(F0, F0[:,1])
print(Fa)
T = 10
F = np.c_[F0, np.zeros((np.shape(F0)[0],T - np.shape(F0)[1]))]
print(np.zeros((3,1)))

A = np.array([[2,3],
              [1,4]])
B = np.array([[1,4],
              [1,1]]).T
C = np.array([[1,4],
              [1,1]])

D = B[:,1] - np.dot(A[0,:],B)

E = np.array([[2,3]]).T

F = np.add(A[:,1],B)

G = np.array([[1, 2, 3, 4, 5], 
              [6, 7, 8, 9, 10]])
print(D)
print(np.dot(C,E))
print(G[:,1])
# A operação de baixo ele duplica o resultado da multiplicação e subtraiu elemento a elemento
print(G[:,1] - np.dot(C,E))
H = np.copy(G[:,1])
print(H - np.dot(C,E))
H = np.zeros((2,1))
H[:,0] = G[:,1]
# OU
H[:,0] = np.copy(G[:,1])
# RESOLVIDO O PROBLEMA!
print(H)
print(H-np.dot(C,E))
# Pra multiplicação não há problema fazer direto a retirada de uma coluna
# Agora, para soma e subtração, é necessário declarar uma variável como matriz
# E alocar a coluna na mesma. 
print(np.dot(C[1,:],E))

I = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
print(I[::2,:]) # pega todas as linhas, de duas em duas
print(I[::2,1]) # pega o segundo elemento de cada linha, de duas em duas linhas

print(np.ones(3)) # um vetor de 3 componentes com 1 em cada componente
print(np.ones((2,3)))

x = np.linspace(0,40,num=2000)
xd = len(x)
print(xd)
y = np.zeros((xd))
alpha = 200
step1 = int(0.05*xd)
print(step1)
y[0:step1] = 0.5*(1+(alpha-1)*x[0:step1]/x[step1])

plt.plot(x,y,'r--')
plt.figure(1,figsize=(8,4))
plt.xlim(0,40)
plt.grid(True)
plt.show()

comp_0 = pd.read_excel('sismo_artificial_comp0.xlsx').to_numpy().T
comp_90 = pd.read_excel('sismo_artificial_comp90.xlsx').to_numpy().T
up = pd.read_excel('sismo_artificial_up.xlsx').to_numpy().T

col = np.shape(comp_0)[1]
sismo_artificial = np.zeros((3,col))

sismo_artificial[0,:] = np.copy(comp_0)
sismo_artificial[1,:] = np.copy(comp_90)
sismo_artificial[2,:] = np.copy(up)

print(sismo_artificial)



