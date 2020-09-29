##### SOLVING VIBRATIONS PROBLEMS BY NUMERICAL METHODS #####

import numpy as np
import numpy.linalg as la
import Dynas as dyn
import pandas as pd
import matplotlib.pyplot as plt

K = np.array([[ 8, -4,  0], 
              [-4,  8, -4],
              [ 0, -4,  4]])
M = np.diag([4, 4, 4])

def modal_analysis(K,M):
    # solving the eigenvalue/eigenvector problem
    D = np.dot(la.inv(M),K)
    lambdak, Phi = la.eig(D)
    # sorting the eigenvalues and eigenvectors in ascending order
    index_lambdak = lambdak.argsort()
    lambdak = lambdak[index_lambdak]
    Phi = Phi[:, index_lambdak]
    # computing the natural angular frequencys and natural frequencys
    wk = np.sqrt(np.real(lambdak))
    fk = wk/(2*np.pi)

    return fk, wk, Phi

def finite_diff(F, x0, v0, dt, M, K, C, T):
    """ SOLVING DIFFERENTIAL EQUATIONS BY THE FINITE DIFFERENCE METHOD
    F  = matrix including in each column the load vector along time with step dt
    x0 = initial position column vector
    v0 = initial velocity column vector
    dt = time step (uniform along duration)
    M  = mass matrix
    K  = stiffness matrix
    C  = damping matrix
    T  = total duration of the analysis (not necessarily the same duration of the load)"""

    ### INITIAL PARAMETERS ####

    # defining the number of steps of analysis = Ns
    Ns =  int(T/dt)+1
    # step t0 (initial acceleration)
    ngl = np.shape(F)[0] # captures the number of degrees of freedom

    ### MODELLING THE DISPLACEMENTS ###

    x_before = np.zeros((ngl,1))
    # matrix that indicates the displacements, in each degree of freedom, along the time of 
    # duration of analysis. Each column is a time step
    x = np.zeros((ngl, Ns))
    x[:,0] = x0[:,0]

    ### SOLVING INITIAL STEP ###

    # initial Force F0 is equivalent to the first column of the matrix of load vectors F along time
    aux1 = np.zeros((ngl,1))
    aux1[:,0] = np.copy(F[:,0])
    aux2 = aux1 - np.dot(C,v0) - np.dot(K,x0)
    a0 = np.dot(la.inv(M),aux2)
    # step t-1 (before initial condition)
    x_before = dt*dt*a0/2 - dt*v0 + x0 
    # step t+1 (after initial condition)
    C1 = M / (dt*dt) + C / (2*dt)
    C2 = K - 2*M / (dt*dt)
    C3 = M / (dt*dt) - C / (2*dt)
    aux3 = aux1 - np.dot(C2, x0) - np.dot(C3, x_before)
    x[:,1] = np.dot(la.inv(C1), aux3[:,0])

    ### INTEGRATING ALONG THE DURATION OS ANALYSIS ###

    i = 0
    aux4 = np.zeros((ngl,1))
    aux5 = np.zeros((ngl,1))
    aux6 = np.zeros((ngl,1))
    aux7 = np.zeros((ngl,1))
    for i in range(1,Ns-1):
        aux4[:,0] = np.copy(F[:,i])
        aux5[:,0] = np.copy(x[:,i])
        aux6[:,0] = np.copy(x[:,i-1])
        aux7[:,0] = np.copy(x[:,i+1])
        aux7 =  np.dot(la.inv(C1), aux4 - np.dot(C2,aux5) - np.dot(C3,aux6))
        x[:,i+1] = np.copy(aux7[:,0])
    return x

t = np.linspace(0, 1, num=11)
F = np.sin(t)

fk, wk, Phi = modal_analysis(K,M)
print(wk[1])
C = dyn.damping_matrix(wk[0], wk[2], 0.01, 0.01, M, K) 
print(C)

F0_teste = np.zeros((2,11))
x0_teste = np.array([[0.5, 0]]).T
v0_teste = np.zeros((2,1))
K_teste = np.array([[ 90,-30],
                    [-30, 30]])
M_teste = np.array([[8,0],
                    [0,4]])
C_teste = np.zeros((2,2))

solution = finite_diff(F0_teste, x0_teste, v0_teste, 0.1, M_teste, K_teste, C_teste, 1)
print(solution.T) 

def Newmark(Kr,Mr,Cr,F,u0,v0,t):
    """Kr : Restricted Stiffness Matrix
       Mr : Restricted Mass Matrix
       Cr : Restricted Damping Matrix
       F  : Load Vector discretized in time
       u0 : initial position vector
       v0 : initial velocity vector
       t  : list with time discretized"""
    tf = int(len(t))
    n = len(F[:,0])
    a = np.zeros((n,tf))
    v = np.zeros((n,tf))
    u = np.zeros((n,tf))
    u[:,0] = np.copy(u0[:,0])
    v[:,0] = np.copy(v0[:,0])
    dt =t[1]-t[0]

    ## DEFINING PARAMETERS OF THE METHOD ##

    delta = 0.5
    alfa = 0.25
    a0 = 1/(alfa*(dt**2))
    a1 = 1/(alfa*dt)
    a2 = (1/(2*alfa))-1
    a3 = (1-delta)*dt
    a4 = delta*dt
    a5 = a4*a0
    a6 = a4*a1-1
    a7 = a4*a2 - a3

    ## WITH EQUATION OF MOVEMENT, u0 AND v0 OBTAIN THE INITIAL ACCELERATION a[:,0] ##
    a[:,0] = np.dot(la.inv(Mr),(F[:,0]-np.dot(Cr,v[:,0])-np.dot(Kr,u[:,0])))
    ## CONSTANT COEFFICIENT IN THE NUMERICAL METHOD ##
    C1 = la.inv((a0*Mr + a5*Cr + Kr))

    for i in range(0,tf-1):
        aux1 = np.dot(Mr,(a0*u[:,i]+ a1*v[:,i] + a2*a[:,i]))
        aux2 = np.dot(Cr,(a5*u[:,i]+ a6*v[:,i] + a7*a[:,i]))
        u_aux = F[:,i+1]+ aux1 + aux2
        u[:,i+1] = np.dot(C1,u_aux)
        v[:,i+1] = a5*(u[:,i+1] - u[:,i]) - a6*v[:,i] - a7*a[:,i]
        a[:,i+1] = a0*(u[:,i+1] - u[:,i]) - a1*v[:,i] - a2*a[:,i]
    
    return u,v,a

desloc, veloc, accel = Newmark(K_teste,M_teste,C_teste,F0_teste,x0_teste,v0_teste,t)
print(desloc.T)
error = desloc.T - solution.T
print("erro relativo aos dois métodos, em % \n", error*100)

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

t = np.linspace(0,39.98,num=2000)

print(t)

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
plt.plot(t,0.01*up,'-b')
plt.title("Componente UP")
plt.xlabel('tempo (s)')
plt.ylabel('Aceleração (m/s²)')
plt.xlim(0,39.98);plt.ylim(-0.3,0.3)
plt.grid(True)


plt.show()
