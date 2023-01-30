import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

plt.style.use('grayscale')

# Pour les intégrations

temps=10000 #s                Intervalle de temps évaluée
pas = 0.01 #s                 Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Données

A=2.01*10**-4 # m^2           Aire de la section transversale
Hp=16.5*10**-3 # m            Hauteur des puits
ho=9.96*10**-3 # m            Hauteur initiale de la colonne de liquide dans les puits
R=10*10**13 # Pa*s/m^3        Résistance hydraulique (considérée constante)
Q1=10*(10**-9)/60 #m^3/s      Débit min
Q2=25*(10**-9)/60 #m^3/s      Débit max 

ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique

#

t_eval = np.arange(0, temps, pas)

# t90 en fonction de ho

# Vecteurs pour les valeurs du débit et du temps t99 calculé pour chaque valeur de ho

hov=[]
t90a=[]
t90b=[]
t90c=[]
Pe=[]

# Boucle pour calculer les temps t90 en faisant varier le débit

for ho in np.arange(0.2, 10, 0.2): #débit en microL/min

    
    # Définition des EDO's à être résolues - c'est juste ho qui change

    Fa = lambda t, H1: Q1/A-Patm/(R*A)*(1-ho*10**-3/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*ho*10**-3/Hp))-2*ro*g/(R*A)*(H1-ho*10**-3)
    Fb = lambda t, H1: Q2/A-Patm/(R*A)*(1-ho*10**-3/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*ho*10**-3/Hp))-2*ro*g/(R*A)*(H1-ho*10**-3)

    
    # Pour chaque valeur de ho, on définit deltaP et l'équation eq qui doit être égale à zéro pour t=t90
    # eq est donc égale au debit Q' qui traverse la biopuce moins 99% du débit Q imposé par la pompe
    
    # Pour ho_a 
    
    def deltaPa(t, H1) :
        return Patm*(1-ho*10**-3/Hp)*(1/(1-H1[0]/Hp)-1/(1+H1[0]/Hp-2*ho*10**-3/Hp))+2*ro*g*(H1[0]-ho*10**-3)
    
    def eqa(t, H1) : 
        return deltaPa(t,H1)/R-0.90*Q1
    
    # Pour ho_b 
    
    def deltaPb(t, H1) :
        return Patm*(1-ho*10**-3/Hp)*(1/(1-H1[0]/Hp)-1/(1+H1[0]/Hp-2*ho*10**-3/Hp))+2*ro*g*(H1[0]-ho*10**-3)
    
    def eqb(t, H1) : 
        return deltaPb(t,H1)/R-0.90*Q2
    
    
    # On résout les EDO's et on trouve les valeurs de t pour lesquels eq=0 avec le paramètre 'events' de la fonction solve_ivp
    
    solsa = solve_ivp(Fa, [0, temps], [ho*10**-3], method='Radau', t_eval=t_eval, events=eqa)
    solsb = solve_ivp(Fb, [0, temps], [ho*10**-3], method='Radau', t_eval=t_eval, events=eqb)

    # A la fin de chaque boucle, on ajoute la valeur de Q et de t99 pour chaque ho dans les vecteurs respectifs 
      
    hov.append(ho)
    t90a.append(solsa.t_events[0]/60)
    t90b.append(solsb.t_events[0]/60)

    

# Courbes de t90 selon Q pour les trois valeurs de ho


plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,10)
#plt.ylim(0,)
plt.plot(hov, t90a, 'o', label='10 microL/min')
plt.plot(hov, t90b, 'o', label='25 microL/min')
plt.xlabel('ho (mm)')
plt.ylabel('t90 (min)')
plt.legend(loc = 7)
plt.show()
