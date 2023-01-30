import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# Pour les intégrations

temps=50000 #s                Intervalle de temps évalué
pas = 1 #s                    Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Données

A=2.01*10**-4 # m^2            Aire de la section transversale
Hp=16.5*10**-3 # m             Hauteur des puits
R=10*10**13 # Pa*s/m^3         Résistance hydraulique (considérée constante)

ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique
Pmax=50*100+Patm # Pa         Pression maximale dans le puits d'entrée

# Les valeurs de ho
    
hoa=4.98*10**-3 #m
hob=7.47*10**-3 #m
hoc=9.96*10**-3 #m

# t90 en fonction du débit Q

# Vecteurs pour les valeurs du débit et du temps t99 calculé pour chaque valeur de ho

Debit=[]
t90a=[]
t90b=[]
t90c=[]
Pe=[]

# Boucle pour calculer les temps t90 en faisant varier le débit

for Q in range (10, 50, 1): #débit en microL/min
    
    Debit.append(Q)

   
    t_eval = np.arange(0, temps, pas)
    
    # Définition et résolution des EDO's pour trouver les instants pour lesquels Pe=Pmax

    Fa = lambda t, H1: Q*(10**-9)/60/A-Patm/(R*A)*(1-hoa/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*hoa/Hp))-2*ro*g/(R*A)*(H1-hoa)
    Fb = lambda t, H1: Q*(10**-9)/60/A-Patm/(R*A)*(1-hob/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*hob/Hp))-2*ro*g/(R*A)*(H1-hob)
    Fc = lambda t, H1: Q*(10**-9)/60/A-Patm/(R*A)*(1-hoc/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*hoc/Hp))-2*ro*g/(R*A)*(H1-hoc)
    

    
    solsa = solve_ivp(Fa, [0, temps], [hoa], method='Radau', t_eval=t_eval)
    solsb = solve_ivp(Fb, [0, temps], [hob], method='Radau', t_eval=t_eval)
    solsc = solve_ivp(Fc, [0, temps], [hoc], method='Radau', t_eval=t_eval)
      
    
    # Boucles pour trouver les instants où Pe=Pmax
    
    for j in range (0, n1, 1) :
        
        if Patm*(Hp-hoa)/(Hp-(solsa.y[0,j])) < Pmax:
            j=j+1
            
        else :
            break
        
    for k in range (0, n1, 1) :
        
        if Patm*(Hp-hob)/(Hp-(solsb.y[0,k])) < Pmax:
            k=k+1
            
        else :
            break
    
    for l in range (0, n1, 1) :
        
        if Patm*(Hp-hoc)/(Hp-(solsc.y[0,l])) < Pmax:
            l=l+1
            
        else :
            break
        
    # Deuxième régime 
    
    Fa1 = lambda t, H: Q/A-Pmax/(R*A)-2*ro*g/(R*A)*(H-hoa)+Patm/(R*A)*(Hp-hoa)/(Hp+H-2*hoa)
    Fb1 = lambda t, H: Q/A-Pmax/(R*A)-2*ro*g/(R*A)*(H-hob)+Patm/(R*A)*(Hp-hob)/(Hp+H-2*hob)
    Fc1 = lambda t, H: Q/A-Pmax/(R*A)-2*ro*g/(R*A)*(H-hoc)+Patm/(R*A)*(Hp-hoc)/(Hp+H-2*hoc)
    
    t_evala1 = np.arange(j*pas, temps, pas)
    t_evalb1 = np.arange(k*pas, temps, pas)
    t_evalc1 = np.arange(l*pas, temps, pas)
    
    
    # Boucles pour trouver les instants où Q'=90%*Q
    
   
    def eqa(t, H1) : 
        return (Pmax-Patm*(Hp-hoa)/(Hp+H1[0]-2*hoa)+2*ro*g*(H1[0]-hoa))/R-0.90*Q*(10**-9)/60

    def eqb(t, H1) : 
        return (Pmax-Patm*(Hp-hob)/(Hp+H1[0]-2*hob)+2*ro*g*(H1[0]-hob))/R-0.90*Q*(10**-9)/60

    def eqc(t, H1) : 
        return (Pmax-Patm*(Hp-hoc)/(Hp+H1[0]-2*hoc)+2*ro*g*(H1[0]-hoc))/R-0.90*Q*(10**-9)/60
    

    solsa1 = solve_ivp(Fa1, [j*pas, temps+1], [hoa], method='Radau', t_eval=t_evala1, events=eqa)
    solsb1 = solve_ivp(Fb1, [k*pas, temps+1], [hob], method='Radau', t_eval=t_evalb1, events=eqb)
    solsc1 = solve_ivp(Fc1, [l*pas, temps+1], [hoc], method='Radau', t_eval=t_evalc1, events=eqc)
    
    t90a.append(j*pas/60+solsa1.t_events[0]/60)
        
    t90b.append(k*pas/60+solsb1.t_events[0]/60)
    
    t90c.append(l*pas/60+solsc1.t_events[0]/60)
        

# Courbes de t99 selon Q pour les trois valeurs de ho


plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(9, 50)
#plt.ylim(0,)
plt.plot(Debit, t90a, 's', color='0.3', label='1,0 mL')
plt.plot(Debit, t90b, 'o', color='0.3', label='1,5 mL')
plt.plot(Debit, t90c, 'd', color='0.3', label='2,0 mL')
plt.xlabel('Q (microL/min)')
plt.ylabel('t90 (min)')
plt.legend(loc = 7, title='Volume initial')
plt.show()
