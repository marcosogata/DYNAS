##### SOLVING VIBRATIONS PROBLEMS BY FINITE DIFFERENCE METHOD #####

import numpy as np
import numpy.linalg as la
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

def damping_matrix(wki, wkj, zti, ztj, M, K):
    # given two pairs (zti, wki) and (ztj, wkj) the alpha coefficients are obtained by solving
    # the linear system between these two pairs
    alpha = la.solve([[1/(2*wki), wki/2], [1/(2*wkj), wkj/2]], [zti, ztj])
    # define damping matrix as a linear combination of K and M
    C = alpha[0]*M + alpha[1]*K
    return C

def finite_diff(F0, x0, v0, dt, M, K, T):
    """ Finite Difference Method
    F0 = matrix including in each column the load vector along time with step dt
    x0 = initial position column vector
    v0 = initial velocity column vector
    dt = time step (uniform along duration)
    M = mass matrix
    K = stiffness matrix
    C = damping matrix
    T = total duration of the analysis (not necessarily the same duration of the load)
    OBS: in the code, the duration of load must be less than the duration of analysis"""

    ### MODELLING THE LOAD ###
    # step t0 (initial acceleration)
    ngl = np.shape(F0)[0] # captures the number of degrees of freedom
    # adding the duration when the force is not being applied anymore
    F = np.c_[F0, np.zeros((np.shape(F0)[0],T - np.shape(F0)[1]))] 

    ### MODELLING THE DISPLACEMENTS ###
    x_before = np.zeros((ngl,1))
    # matrix that indicates the displacements, in each degree of freedom, along the time of 
    # duration of analysis. Each column is a time step
    x = np.zeros((ngl, T))

    ### SOLVING INITIAL STEP ###
    # initial Force F0 is equivalent to the first column of the matrix of load vectors F along time
    Ma0 = F[:,0] - np.dot(C,v0) - np.dot(K,x0)
    a0 = np.dot(la.inv(M), Ma0)
    # step t-1 (before initial condition)
    x_before = dt*dt*a0/2 - dt*v0 + x0 
    # step t+1 (after initial condition)
    C1 = M / (dt*dt) + C / (2*dt)
    C2 = K - 2*M / (dt*dt)
    C3 = M / (dt*dt) - C / (2*dt)
    x[:,1] = np.dot(la.inv(C1), F0 - np.dot(C2, x0) - np.dot(C3, x_before))

    ### INTEGRATING ALONG THE DURATION OS ANALYSIS ###
    i = 0
    for i in range(0,T+1):
        i = i+1
        x[:,i+1] = np.dot(la.inv(C1), F[:,i] - np.dot(C2, x[:,i]) - np.dot(C3, x[:,i-1]))

    return x

t = np.linspace(0,10,1)
F = np.sin(t)





    pass

D = modal_analysis(K,M)
C = damping_matrix(D[1][0], D[1][2], 0.01, 0.01, M, K) 

print(C)

