import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

plt.style.use('grayscale')

# Pour les intégrations

temps=3000 #s                 Intervalle de temps évaluée
pas = 0.01 #s                 Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Données

A=2.01*10**-4 # m^2           Aire de la section transversale
Hp=16.5*10**-3 # m            Hauteur des puits
ho=9.96*10**-3 # m            Hauteur initiale de la colonne de liquide dans les puits
R=10*10**13 # Pa*s/m^3        Résistance hydraulique (considérée constante)
Q2=10*(10**-9)/60 #m^3/s      Débit max 

ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique

# Les valeurs de ho

hoa=4.98*10**-3 #m
hob=7.47*10**-3 #m
hoc=9.96*10**-3 #m

t_eval = np.arange(0, temps, pas)

# t90 en fonction du débit Q

# Vecteurs pour les valeurs du débit et du temps t99 calculé pour chaque valeur de ho

Debit=[]
t90a=[]
t90b=[]
t90c=[]
Pe=[]

# Boucle pour calculer les temps t90 en faisant varier le débit

for Q in range (1, 50, 1): #débit en microL/min

    
    # Définition des EDO's à être résolues - c'est juste ho qui change

    Fa = lambda t, H1: Q*(10**-9)/60/A-Patm/(R*A)*(1-hoa/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*hoa/Hp))-2*ro*g/(R*A)*(H1-hoa)
    Fb = lambda t, H1: Q*(10**-9)/60/A-Patm/(R*A)*(1-hob/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*hob/Hp))-2*ro*g/(R*A)*(H1-hob)
    Fc = lambda t, H1: Q*(10**-9)/60/A-Patm/(R*A)*(1-hoc/Hp)*(1/(1-H1/Hp)-1/(1+H1/Hp-2*hoc/Hp))-2*ro*g/(R*A)*(H1-hoc)
    
    # Pour chaque valeur de ho, on définit deltaP et l'équation eq qui doit être égale à zéro pour t=t90
    # eq est donc égale au debit Q' qui traverse la biopuce moins 99% du débit Q imposé par la pompe
    
    # Pour ho_a 
    
    def deltaPa(t, H1) :
        return Patm*(1-hoa/Hp)*(1/(1-H1[0]/Hp)-1/(1+H1[0]/Hp-2*hoa/Hp))+2*ro*g*(H1[0]-hoa)
    
    def eqa(t, H1) : 
        return deltaPa(t,H1)/R-0.90*Q*(10**-9)/60
    
    # Pour ho_b 
    
    def deltaPb(t, H1) :
        return Patm*(1-hob/Hp)*(1/(1-H1[0]/Hp)-1/(1+H1[0]/Hp-2*hob/Hp))+2*ro*g*(H1[0]-hob)
    
    def eqb(t, H1) : 
        return deltaPb(t,H1)/R-0.90*Q*(10**-9)/60
    
    # Pour ho_c 
    
    def deltaPc(t, H1) :
        return Patm*(1-hoc/Hp)*(1/(1-H1[0]/Hp)-1/(1+H1[0]/Hp-2*hoc/Hp))+2*ro*g*(H1[0]-hoc)
    
    def eqc(t, H1) : 
        return deltaPc(t,H1)/R-0.90*Q*(10**-9)/60
    
    # On résout les EDO's et on trouve les valeurs de t pour lesquels eq=0 avec le paramètre 'events' de la fonction solve_ivp
    
    solsa = solve_ivp(Fa, [0, temps], [hoa], method='Radau', t_eval=t_eval, events=eqa)
    solsb = solve_ivp(Fb, [0, temps], [hob], method='Radau', t_eval=t_eval, events=eqb)
    solsc = solve_ivp(Fc, [0, temps], [hoc], method='Radau', t_eval=t_eval, events=eqc)

    # A la fin de chaque boucle, on ajoute la valeur de Q et de t99 pour chaque ho dans les vecteurs respectifs 
      
    Debit.append(Q)
    t90a.append(solsa.t_events[0]/60)
    t90b.append(solsb.t_events[0]/60)
    t90c.append(solsc.t_events[0]/60)
    
    # On calcule aussi la valeur de la pression relative dans le puits d'entrée pour chaque valeur de débit évaluée
    
    def Pe_a(t, H1) :
        return Patm*(Hp-hoa)/(Hp-H1)-Patm
    
    Pe.append(Pe_a(solsa.t_events[0], solsa.y[0,n1])*10**-2)

# Courbes de t90 selon Q pour les trois valeurs de ho


plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,50)
#plt.ylim(0,)
plt.plot(Debit, t90a, 's', color='0.3', label='1,0 mL')
plt.plot(Debit, t90b, 'o', color='0.3', label='1,5 mL')
plt.plot(Debit, t90c, 'd', color='0.3', label='2,0 mL')
plt.xlabel('Q (microL/min)')
plt.ylabel('t90 (min)')
plt.legend(loc = 7)
plt.show()

# Courbe de Pe selon Q 

plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,50)
plt.ylim(0, 500)
plt.plot(Debit, Pe, 'o')
plt.title("Pression relative dans le puits d'entrée en t99 en fonction du débit")
plt.xlabel('Q (microL/min)')
plt.ylabel("Pe relative (mbar)")
plt.show()

