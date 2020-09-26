import numpy as np
import numpy.linalg as la
import scipy.linalg as sc
import matplotlib.pyplot as plt
import pandas as pd
from pandas import ExcelWriter
from numba import jit
from mpl_toolkits.mplot3d import Axes3D 

##############################################
### BUILDING THE STIFFNESS AND MASS MATRIX ###
##############################################

def Matrix3D(Name,Nn):
    """ METHOD FOR BUILDING STIFFNESS AND MASS MATRIX FOR 3D FRAME STRUCTURE
        Name : File Name (remember to put ".xlsx").
        Nn : Number of nodes.
        OBS: this method only works if the worksheet is exactly formatted as the file FRAME3D.xlsx"""

    File = pd.read_excel(Name)

    ### BUILDING THE COORDINATE ROW VECTORS (X AND Y FOR EACH NODE) ###

    cx = list(File['Cx'])[0:Nn]
    cy = list(File['Cy'])[0:Nn]
    cz = list(File['Cz'])[0:Nn]
    
    ### BUILDING THE MATRIX ID IN NODAL PARTITION (INITIAL AND END NODE OF EACH BAR) ###
    ### Matrix ID - identify the connection between nodes         ###

    # connectivities of each bar
    Id1 = list(File['barra (nó 1)'])
    Id2 = list(File['barra (nó 2)'])
    # number of bars
    Nb = len(Id1) 
    # Build a matrix 2xNb were the first row corresponds to the i node and the second to the j node
    ID =np.zeros((2,Nb))
    # allocating the i node and j node
    ID[0,:] = Id1
    ID[1,:] = Id2

    ### ALLOCATE THE PROPERTIES OF THE SECTIONS IN EACH BAR ###

    A = list(File['Area(m2)'])
    Iz =list(File['Inércia z (m4)'])
    Iy =list(File['Inércia y (m4)'])
    RHO = list(File['Densidade'])
    
    ### PROPERTIES OF THE MATERIAL USING SI UNITS ###

    E = 28*10**9 
    G = E / 2.4
    J = Iy + Iz

    ###  BUILDING THE MATRIX ID IN TERMS OF DEGREES OF FREEDOM ###
    ###  Each node has six degrees of freedom                  ###

    # first row of matrix ID desmembered in the six degrees of freedom per node
    lb = ID[0,:]
    l1 = lb*6-5
    l2 =lb*6-4
    l3 =lb*6-3
    l4 = lb*6-2
    l5 =lb*6-1
    l6 =lb*6
    # second row of matrix ID desmembered in the six degrees of freedom per node
    lb2 =ID[1,:]
    l7 = lb2*6-5
    l8 =lb2*6-4
    l9 =lb2*6-3
    l10 = lb2*6-2
    l11 =lb2*6-1
    l12 =lb2*6
    # build the (12 degrees of freedom per element)x(Number of bars) matrix
    IDG = np.zeros((12,Nb))
    # allocating each degree of freedom in each row of the matrix ID
    IDG[0,:]=l1
    IDG[1,:]=l2
    IDG[2,:]=l3
    IDG[3,:]=l4
    IDG[4,:]=l5
    IDG[5,:]=l6
    IDG[6,:]=l7
    IDG[7,:]=l8
    IDG[8,:]=l9
    IDG[9,:]=l10
    IDG[10,:]=l11
    IDG[11,:]=l12

    ### BUILDING THE STIFFNESS AND MASS MATRIX OF THE 3D BAR ELEMENT ###

    # Size of the matrix is equal to the number of nodes multiply by the number of degrees of freedom per node
    K = np.zeros((Nn*6,Nn*6))
    M = np.zeros((Nn*6,Nn*6))

    for i in range (Nb):
        ## CAPTURING THE NODES BASED ON THE MATRIX ID IN NODAL PARTITION
        k1 = int(ID[0,i] -1)
        k2 = int(ID[1,i] -1)
        ## CALCULATING THE LENGTH OF EACH ELEMENT 
        Lx = cx[k2] - cx[k1]
        Ly = cy[k2] - cy[k1]
        Lz = cz[k2] - cz[k1]
        L = np.sqrt(Lx**2 + Ly**2 +Lz**2)
        ## CALCULATING THE DIRECTION COSINES
        l = Lx/L
        m = Ly/L
        n = Lz/L
        ## DEFINING THE ROTATIONAL MATRIX (ONLY FOR ORTHOGONAL FRAMES)
        Cxz = np.sqrt(l**2 + n**2)
        if Cxz ==0:
            R = np.array([[ 0,m,0, 0,0,0, 0,0,0, 0,0,0],
                          [-m,0,0, 0,0,0, 0,0,0, 0,0,0],
                          [ 0,0,1, 0,0,0, 0,0,0, 0,0,0],
                          [ 0,0,0, 0,m,0, 0,0,0, 0,0,0],
                          [ 0,0,0,-m,0,0, 0,0,0, 0,0,0],
                          [ 0,0,0, 0,0,1, 0,0,0, 0,0,0],
                          [ 0,0,0, 0,0,0, 0,m,0, 0,0,0],
                          [ 0,0,0, 0,0,0,-m,0,0, 0,0,0],
                          [ 0,0,0, 0,0,0, 0,0,1, 0,0,0],
                          [ 0,0,0, 0,0,0, 0,0,0, 0,m,0],
                          [ 0,0,0, 0,0,0, 0,0,0,-m,0,0],
                          [ 0,0,0, 0,0,0, 0,0,0, 0,0,1]])
        else:
            R = np.array([[       l,    m,        n,         0,    0,        0,         0,     0,        0,         0,    0,        0],
                          [-l*m/Cxz,  Cxz, -m*n/Cxz,         0,    0,        0,         0,     0,        0,         0,    0,        0],
                          [  -n/Cxz,    0,    l/Cxz,         0,    0,        0,         0,     0,        0,         0,    0,        0],
                          [       0,    0,        0,         l,    m,        n,         0,     0,        0,         0,    0,        0],
                          [       0,    0,        0,  -l*m/Cxz,  Cxz, -m*n/Cxz,         0,     0,        0,         0,    0,        0],
                          [       0,    0,        0,    -n/Cxz,    0,    l/Cxz,         0,     0,        0,         0,    0,        0],
                          [       0,    0,        0,         0,    0,        0,         l,     m,        n,         0,    0,        0],
                          [       0,    0,        0,         0,    0,        0,  -l*m/Cxz,   Cxz, -m*n/Cxz,         0,    0,        0],
                          [       0,    0,        0,         0,    0,        0,    -n/Cxz,     0,    l/Cxz,         0,    0,        0],
                          [       0,    0,        0,         0,    0,        0,         0,     0,        0,         l,    m,        n],
                          [       0,    0,        0,         0,    0,        0,         0,     0,        0,  -l*m/Cxz,  Cxz, -m*n/Cxz],
                          [       0,    0,        0,         0,    0,        0,         0,     0,        0,    -n/Cxz,    0,    l/Cxz]]) 
   
    ### LOCAL STIFFNESS MATRIX OF THE ELEMENT ###

        Ke =np.array(([[ E*A[i]/L,                 0,                 0,        0,                0,                0,-E*A[i]/L,                 0,                 0,        0,                0,                0],
                       [        0, 12*E*Iz[i]/(L**3),                 0,        0,                0, 6*E*Iz[i]/(L**2),        0,-12*E*Iz[i]/(L**3),                 0,        0,                0, 6*E*Iz[i]/(L**2)],
                       [        0,                 0, 12*E*Iy[i]/(L**3),        0,-6*E*Iy[i]/(L**2),                0,        0,                 0,-12*E*Iy[i]/(L**3),        0,-6*E*Iy[i]/(L**2),                0],
                       [        0,                 0,                 0, G*J[i]/L,                0,                0,        0,                 0,                 0,-G*J[i]/L,                0,                0],
                       [        0,                 0, -6*E*Iy[i]/(L**2),        0,      4*E*Iy[i]/L,                0,        0,                 0,  6*E*Iy[i]/(L**2),        0,      2*E*Iy[i]/L,                0],
                       [        0,  6*E*Iz[i]/(L**2),                 0,        0,                0,      4*E*Iz[i]/L,        0, -6*E*Iz[i]/(L**2),                 0,        0,                0,      2*E*Iz[i]/L],
                       [-E*A[i]/L,                 0,                 0,        0,                0,                0, E*A[i]/L,                 0,                 0,        0,                0,                0],
                       [        0,-12*E*Iz[i]/(L**3),                 0,        0,                0,-6*E*Iz[i]/(L**2),        0, 12*E*Iz[i]/(L**3),                 0,        0,                0,-6*E*Iz[i]/(L**2)],
                       [        0,                 0,-12*E*Iy[i]/(L**3),        0, 6*E*Iy[i]/(L**2),                0,        0,                 0, 12*E*Iy[i]/(L**3),        0, 6*E*Iy[i]/(L**2),                0],
                       [        0,                 0,                 0,-G*J[i]/L,                0,                0,        0,                 0,                 0, G*J[i]/L,                0,                0],
                       [        0,                 0, -6*E*Iy[i]/(L**2),        0,      2*E*Iy[i]/L,                0,        0,                 0,  6*E*Iy[i]/(L**2),        0,      4*E*Iy[i]/L,                0],
                       [        0,  6*E*Iz[i]/(L**2),                 0,        0,                0,      2*E*Iz[i]/L,        0, -6*E*Iz[i]/(L**2),                 0,        0,                0,      4*E*Iz[i]/L]]))
    
    ### LOCAL MASS MATRIX OF THE ELEMENT ###
    
        Jx = J[i]/A[i]
        Me = ((RHO[i]*A[i]*L)/420)*np.array([[140,    0,    0,      0,      0,      0,  70,    0,    0,      0,      0,      0],
                                             [  0,   156,   0,      0,      0,   22*L,   0,   54,    0,      0,      0,  -13*L],
                                             [  0,    0,  156,      0,  -22*L,      0,   0,    0,   54,      0,   13*L,      0],
                                             [  0,    0,    0, 140*Jx,      0,      0,   0,    0,    0,  70*Jx,      0,      0],
                                             [  0,    0,-22*L,      0, 4*L**2,      0,   0,    0,-13*L,      0,-3*L**2,      0],
                                             [  0, 22*L,    0,      0,      0, 4*L**2,   0, 13*L,    0,      0,      0,-3*L**2],
                                             [ 70,    0,    0,      0,      0,      0, 140,    0,    0,      0,      0,      0],
                                             [  0,   54,    0,      0,      0,   13*L,   0,  156,    0,      0,      0,  -22*L],
                                             [  0,    0,   54,      0,  -13*L,      0,   0,    0,  156,      0,   22*L,      0],
                                             [  0,    0,    0,  70*Jx,      0,      0,   0,    0,    0, 140*Jx,      0,      0],
                                             [  0,    0, 13*L,      0,-3*L**2,      0,   0,    0, 22*L,      0, 4*L**2,      0],
                                             [  0,-13*L,    0,      0,      0,-3*L**2,   0,-22*L,    0,      0,      0, 4*L**2]])
   
   
        ### ROTATION OF THE MATRICES TO THE GLOBAL REFFERENCE SYSTEM ###    
                     
        KT = np.dot(np.dot(R.T, Ke),R) 
        MT = np.dot(np.dot(R.T, Me),R)
        
        ### TEMPORARY MATRICES FOR ALLOCATION IN THE GLOBAL MATRIX OF THE STRUCTURE ###

        k_temp1 = np.zeros((Nn*6,Nn*6))
        m_temp1 = np.zeros((Nn*6,Nn*6))

        ### ALLOCATION OF THE TEMPORARY MATRICES IN THE GLOBAL MATRIX ###

        # Identifying the index os degrees of freedom, to be allocated in the global matrix
        j= int(IDG[0,i]-1)
        f= int(IDG[5,i])
        o= int(IDG[6,i]-1)
        p= int(IDG[11,i])
    
        # constructing the temporary matrices to allocate in the global matrix
        k_temp1[j:f,j:f] = KT[0:6,0:6]
        k_temp1[o:p,j:f] = KT[6:12,0:6]
        k_temp1[j:f,o:p] = KT[0:6,6:12]
        k_temp1[o:p,o:p] = KT[6:12,6:12]
        #Summation including the local matrix in the global matrix
        K += k_temp1 

        #Same procedure applied to the mass matrix
        m_temp1[j:f,j:f] = MT[0:6,0:6]
        m_temp1[o:p,j:f] = MT[6:12,0:6]
        m_temp1[j:f,o:p] = MT[0:6,6:12]
        m_temp1[o:p,o:p] = MT[6:12,6:12]
        M += m_temp1 
        
    return K,M

def Restr(K, M, Nr): 
    """ M  : Mass matrix
        K  : Stiffness matrix
        Nr : List with restricted degrees of freedom
        Obs: The damping matrix considered is Rayleigh Damping Matrix, which consists in
        a linear combination of M and K. Hence, C is consequence of Mr and Kr"""

    ### EXCLUDING THE RESTRICTED ROWS AND COLUMNS

    Kr_1 = np.delete(K,Nr,0)
    Kr   = np.delete(Kr_1,Nr,1)
    Mr_1 = np.delete(M,Nr,0)
    Mr   = np.delete(Mr_1,Nr,1)

    ### EXPORT THE MATRICES TO AN EXCEL FILE

    df1    = pd.DataFrame(Mr)
    writer = ExcelWriter('Matriz de Massa.xlsx')
    df1.to_excel(writer,'Sheet1', index=False)
    writer.save()

    df     = pd.DataFrame(Kr)
    writer = ExcelWriter('Matriz de Rigidez.xlsx')
    df.to_excel(writer,'Sheet1', index=False) 
    writer.save()
 
    return Kr, Mr


def damping_matrix(wki, wkj, zti, ztj, Mr, Kr):
    # given two pairs (zti, wki) and (ztj, wkj) the alpha coefficients are obtained by solving
    # the linear system between these two pairs
    alpha = la.solve([[1/(2*wki), wki/2], [1/(2*wkj), wkj/2]], [zti, ztj])
    # define damping matrix as a linear combination of K and M
    Cr = alpha[0]*Mr + alpha[1]*Kr

    ### EXPORTING TO EXCEL ###

    df     = pd.DataFrame(Cr)
    writer = ExcelWriter('Matriz de Amortecimento.xlsx')
    df.to_excel(writer,'Sheet1', index=False) 
    writer.save()

    return Cr

##############################################
###             DEFINING LOADS             ###
##############################################

### ESPECTRO DE KANAI-TAJIN ###

def Kanai_Tajimi(Ap,tipo,duraçao,dt):
    ## Ap : Peek Ground Acceleration
    ## tipo : tipo de solo 
    ## duraçao : tempo de duração do sinal
    ## dt : Passo temporal
    g = 980.6
    Ap *= g
    pg = 3
    tf = int(duraçao/dt)
    if tipo == '90 DEG': #CONSIDERA UM PGA= -65.84cm/sec/sec
        wg = 16.03 #rad/s
        zg = 0.56
    elif tipo == 'UP':#CONSIDERA UM PGA= 26.95cm/sec/sec
        wg = 16.19 #rad/s
        zg = 0.56
    elif tipo == '0 DEG':#CONSIDERA UM PGA= 28.08cm/sec/sec
        wg = 16.07 #rad/s
        zg = 0.58
        
    f = np.linspace(0,25,tf)
    df = f[1]-f[0]
    w = 2*np.pi*f
    S0 = (Ap**2)/((pg**2)*(np.pi*wg*((1/(2*zg))+2*zg)))
    Sg = S0*((1+4*(zg**2)*(w/wg)**2)/(((1-(w/wg)**2)**2)+4*(zg**2)*(w/wg)**2))


    plt.figure(2, figsize=(8,4)) 
    plt.plot(f,Sg,'b')
    plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)');
    plt.xlim(0,20); plt.ylim(0,max(Sg)*2);plt.title(' Espectro de aceleração')
    plt.grid(True)

    import random 
    P = np.zeros(tf)
    for i in range(tf):
        P1 = random.uniform(0,2*np.pi)
        P[i] = P1
   
    ag= np.zeros(tf)
    t = np.linspace(0,duraçao,tf)
    S = np.zeros(tf)



    for i in range(tf):

        S =np.sqrt(2*Sg*df)*np.cos(w*t[i]+ P)
        ag[i] = sum(S)

    ag*= Ap/np.max(abs(ag))  ## Normalização das acelerações
    
## Função de envoltória

    env =np.ones(tf)
    env1 = np.arange(int(0.05*tf))/(0.05*tf)
    env2 = (1.148698355**t[0:int(0.8*tf)])/64/4
    env[0:int(0.05*tf)] = env1
    env[int(0.2*tf):tf] = env2[::-1]

    plt.figure(3,figsize=(8,4))

    plt.plot(t,ag,'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (cm/s²)');
    plt.xlim(0,duraçao); plt.ylim(-100000,100000);plt.title(' Aceleração do solo')
    plt.grid(True)

    plt.figure(4,figsize=(8,4))
    plt.plot(t,ag,'c')
    plt.plot(t,env,'r--',t,-env,'r--')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (cm/s²)');
    plt.xlim(0,duraçao); plt.ylim(-2*Ap,2*Ap);plt.title(' Função de envoltória')
    plt.grid(True)

    age = ag*env
    age*= Ap/np.max(abs(age))
    
    plt.figure(5,figsize=(8,4))
    plt.plot(t,age,'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (cm/s²)');
    plt.xlim(0,duraçao); plt.ylim(-Ap,Ap);
    plt.title('Aceleração do solo parametrizada')
    plt.grid(True)
    
    return t,age




def Seismic3D(Mr,age,t):
    """Mr : mass matrix restricted 
       t  : list with discrete time 
       age: acceleration signal discretized"""
    tf = int(len(t))
    n = int (len(Mr[0,:]))
    B = np.zeros(n)
    B[::6] = np.ones(int(n/6))
    F = np.zeros((n,tf))
    F1 = np.dot(Mr,B)

    for i in range (tf):
        F[:,i] = F1* age[i] 
    
    plt.figure(6,figsize=(8,4))   
    plt.plot(t,F[n-6 ,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Força N')
    plt.xlim(0,max(t))
    plt.ylim(-max(F[n-6])*1.5,max(F[n-6])*1.5)
    plt.title('Força no 10º pavimento')
    plt.grid(True)
    
    return F

##############################################
### MODAL ANALYSIS AND NUMERICAL SOLUTIONS ###
##############################################

### MODAL ANALYSIS ###

def modal_analysis(K,M,N):
    """K  : Stiffness Matrix 
       M  : Mass Matrix     
       N  : Number of modes desired"""

    # solving the eigenvalue/eigenvector problem
    D = np.dot(la.inv(M),K)
    lambdak, Phi = la.eig(D)
    # sorting the eigenvalues and eigenvectors in ascending order
    index_lambdak = lambdak.argsort()
    lambdak = lambdak[index_lambdak]
    Phi = Phi[:, index_lambdak]
    # computing the natural angular frequencys and natural frequencys
    wk = np.sqrt((lambdak))
    fk = wk/(2*np.pi)

    for k in range(N):
        print(k+1,f"ª Frequência Natural = {fk[k]:3.2f}Hz","\n")

    return fk, wk, Phi

### FINITE DIFFERENCE METHOD ###

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

### NEWMARK METHOD FOR NUMERICAL SOLUTIONS ###

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

