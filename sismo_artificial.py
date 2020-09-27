import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### TRATAMENTO DOS DADOS DE SISMO ### 

# PGA 90 DEG = 65.843 cm/s/s
# PGA 0  DEG = 28.081 cm/s/s
# PGA UP     = 26.952 cm/s/s

### ESPECTRO DE KANAI-TAJIN ###

def Kanai_Tajimi(Name,Ap,tipo,duraçao,dt):
    ## Ap : Peek Ground Acceleration, in % of gravity acceleration g
    ## tipo : soil type 
    ## duraçao : time duration of the signal
    ## dt : step time
    g = 980.6
    Ap *= g
    pg = 3
    tf = int(duraçao/dt)

    ### OBTAINING THE GROUND PROPERTIES WITH RESPECT TO PGA ###

    if tipo == '90 DEG': # for a PGA= -65.843 cm/sec/sec
        wg = 16.03 # rad/s
        zg = 0.56
    elif tipo == 'UP':# for a PGA= 26.952 cm/sec/sec
        wg = 16.19 # rad/s
        zg = 0.56
    elif tipo == '0 DEG':# for a PGA= 28.081 cm/sec/sec
        wg = 16.07 # rad/s
        zg = 0.58
    
    f = np.linspace(0,25,tf)
    df = f[1]-f[0]
    w = 2*np.pi*f
    S0 = (Ap**2)/((pg**2)*(np.pi*wg*((1/(2*zg))+2*zg)))
    Sg = S0*((1+4*(zg**2)*(w/wg)**2)/(((1-(w/wg)**2)**2)+4*(zg**2)*(w/wg)**2))

    plt.figure(2, figsize=(12,4)) 
    plt.plot(f,Sg,'b')
    plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(cm²/s³)');
    plt.xlim(0,20); plt.ylim(0,max(Sg)*1.25);plt.title(' Espectro de aceleração')
    plt.grid(True)

    import random 
    P = np.zeros(tf)
    for i in range(tf):
        P1 = random.uniform(0,2*np.pi)
        P[i] = P1
   
    ag= np.zeros(tf)
    t = np.linspace(0,duraçao,tf)
    S = np.zeros(tf)

    ## OBTAINING THE SIGNAL IN THE TIME DOMAIN (REFERENCE TO SHINOZUKA & JAN (1972))

    for i in range(tf):

        S =np.sqrt(2*Sg*df)*np.cos(w*t[i]+ P)
        ag[i] = sum(S)

    ag*= Ap/np.max(abs(ag))  # Normalization of the accelerations
    
## Envelope function by Barone

    env =np.ones(tf)
    alpha = 10      # arbitrary coefficient
    k1 = -2.5    # arbitrary coefficient
    step1 = int(0.05*tf) # 5% of the total duration
    step2 = int(0.20*tf) # 20% of the total duration
    
    for i in range(0,step1):
        env[i] = S0*(1+(alpha-1)*(t[i]/t[step1]))
    env[step1:step2] = S0*alpha
    for i in range(step2,tf):
        env[i] = S0*alpha*(t[i]/t[step2])**k1

    plt.figure(3,figsize=(12,4))

    plt.plot(t,ag,'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (cm/s²)')
    plt.xlim(0,duraçao); plt.ylim(-1.1*Ap,1.1*Ap);plt.title(' Aceleração do solo')
    plt.grid(True)

    plt.figure(4,figsize=(12,4))
    plt.plot(t,ag,'c')
    plt.plot(t,env,'r--',t,-env,'r--')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (cm/s²)')
    plt.xlim(0,duraçao); plt.ylim(-1.1*Ap,1.1*Ap);plt.title(' Função de envoltória')
    plt.grid(True)

    age = ag*env
    age*= Ap/np.max(abs(age))
    
    plt.figure(5,figsize=(12,4))
    plt.plot(t,age,'g')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (cm/s²)')
    plt.xlim(0,duraçao); plt.ylim(-1.1*Ap,1.1*Ap)
    plt.title('Aceleração do solo parametrizada')
    plt.grid(True)

    plt.show()

    df     = pd.DataFrame(age)
    writer = pd.ExcelWriter(f'{Name}.xlsx')
    df.to_excel(writer,'Sheet1', index=False) 
    writer.save()
    
    return t,age

t90, age90 = Kanai_Tajimi('sismo_artificial_comp90',-65.843/980.6, '90 DEG', 40, 0.02)
t0, age0 = Kanai_Tajimi('sismo_artificial_comp0',28.081/980.6, '0 DEG', 40, 0.02)
tup, ageup = Kanai_Tajimi('sismo_artificial_up',26.952/980.6, 'UP', 40, 0.02)