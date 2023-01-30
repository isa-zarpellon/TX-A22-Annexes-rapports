import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

plt.style.use('grayscale')


# Données

A=2.01*10**-4 # m^2           Aire de la section transversale
ho=9.95*10**-3 # m            Hauteur initiale de la colonne de liquide dans les puits
Hp=16.5*10**-3 # m            Hauteur des puits
R=10*10**13 # Pa*s/m^3        Résistance hydraulique (considérée constante)

Q2=25*(10**-9)/60 # m^3/s     Débit 2
ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique
Pmax=50*100+Patm # Pa         Pression maximale dans le puits d'entrée

# Pour les intégrations

temps=40000 #s                Intervalle de temps évalué
pas = 1 #s                    Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Définition des EDO's 

F2 = lambda t, H: Q2/A-Patm/(R*A)*(1-ho/Hp)*(1/(1-H/Hp)-1/(1+H/Hp-2*ho/Hp))-2*ro*g/(R*A)*(H-ho)

# Résolution des EDO's en utilisant le solver 
# Méthode Runge-Kutta implicite ordre 5


t_eval = np.arange(0, temps, pas)
sol2 = solve_ivp(F2, [0, temps], [ho], method='Radau', t_eval=t_eval, dense_output=True)


# #_________________________________________________________________________________________

# Modèle avec fuite dans le puits d'entrée - résolu pour Q2

# Détermination de l'instant où Pe >= Pmax


for j in range (0, n1, 1) :
    
    if Patm*(Hp-ho)/(Hp-(sol2.y[0,j])) < Pmax:
        j=j+1
        
    else :
        break

# EDO qui définit le deuxième régime

F22 = lambda t, H: Q2/A-Pmax/(R*A)-2*ro*g/(R*A)*(H-ho)+Patm/(R*A)*(Hp-ho)/(Hp+H-2*ho)

# Les nouveaux intervalles de temps

t_evala = np.arange(0, j*pas, pas)
t_evalb = np.arange(j*pas, temps, pas)
     
# Recalcul de la première solution avec le nouvel intervalle 
# Cacul de la deuxième solution

sol2a = solve_ivp(F2, [0, j*pas+1], [ho], method='Radau', t_eval=t_evala)
sol2b = solve_ivp(F22, [j*pas, temps], [sol2.y[0,j]], method='Radau', t_eval=t_evalb)

# DeltaP

deltaP2a = Patm*(1-ho/Hp)*(1/(1-sol2a.y/Hp)-1/(1+sol2a.y/Hp-2*ho/Hp))+2*ro*g*(sol2a.y-ho)
deltaP2b = Pmax-Patm*(Hp-ho)/(Hp+sol2b.y-2*ho)+2*ro*g*(sol2b.y-ho)

# Courbe de H en fonction du temps
Hmax=Hp*np.ones(np.shape(sol2.t))

plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,temps)
#plt.ylim(0,0.6)
plt.plot(sol2.t, (sol2.y[0])*1000, linestyle='solid', color='0', linewidth=1.5, label="Modèle sans fuite")
# plt.plot(sol2.t, (Hmax)*1000, 'y', label="Hmax")
plt.plot(sol2a.t, (sol2a.y[0])*1000, linestyle='dashdot', color='0', linewidth=1.5, label="Modèle avec fuite dans le puits d'entrée")
plt.plot(sol2b.t, (sol2b.y[0])*1000, linestyle='dashdot', color='0', linewidth=1.5)
plt.xlabel('t (s)')
plt.ylabel('H (mm)')
plt.legend(loc = 4)
plt.show()

# Courbe de Ps en fonction du temps


plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,temps)
#plt.ylim(0,0.6)
plt.plot(sol2.t, (Patm*(Hp-ho)/(Hp+sol2.y[0]-2*ho)-Patm)*10**-2, linestyle='solid', color='0', linewidth=1.5, label="Modèle sans fuite" )
plt.plot(sol2a.t, (Patm*(Hp-ho)/(Hp+sol2a.y[0]-2*ho)-Patm)*10**-2, linestyle='dashdot', color='0', linewidth=1.5, label="Modèle avec fuite dans le puits d'entrée")
plt.plot(sol2b.t, (Patm*(Hp-ho)/(Hp+sol2b.y[0]-2*ho)-Patm)*10**-2, linestyle='dashdot', color='0', linewidth=1.5)
plt.ylabel('Ps rel (mbar)')
plt.legend(loc = 1)
plt.show()

# Courbe de Pe en fonction du temps

plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,temps)
#plt.ylim(0,60)
plt.plot(sol2.t, (Patm*(Hp-ho)/(Hp-sol2.y[0])-Patm)*10**-2, linestyle='solid', linewidth=1.5, label="Modèle sans fuite" )
plt.plot(sol2a.t, (Patm*(Hp-ho)/(Hp-sol2a.y[0])-Patm)*10**-2, linestyle='dashdot', color='0', linewidth=1.5, label="Modèle avec fuite dans le puits d'entrée")
plt.plot(sol2b.t, (Pmax-Patm)*(10**-2)*np.ones(np.shape(sol2b.t)), color='0', linewidth=1.5, linestyle='dashdot')
plt.xlabel('t (s)')
plt.ylabel('Pe relative (mbar)')
plt.legend(loc = 4)
plt.show()


#________________________________________________________________________________________
# Calcul de Heq en faisant dH/dt=0

vPmax=[]
vHeq1=[]
vHeq2=[]

R1=5*10**12 # Pa*s/m^3
R2=1*10**14 # Pa*s/m^3

for Pmax1 in np.arange (0, 90, 2) :
   
    a1=(Pmax1*100+Patm-Q2*R1)/(2*ro*g)+Hp-3*ho
    b1=((Q2*R1-(Pmax1*100+Patm))/(2*ro*g)+ho)*(2*ho-Hp)+Patm/(2*ro*g)*(ho-Hp)
    
    a2=(Pmax1*100+Patm-Q2*R2)/(2*ro*g)+Hp-3*ho
    b2=((Q2*R2-(Pmax1*100+Patm))/(2*ro*g)+ho)*(2*ho-Hp)+Patm/(2*ro*g)*(ho-Hp)

    f1 = lambda Hi: Hi**2 + a1*Hi + b1
    f2 = lambda Hi: Hi**2 + a2*Hi + b2

    Heq1 = fsolve(f1, Hp)
    Heq2 = fsolve(f2, Hp)
    
    vPmax.append(Pmax1)
    vHeq1.append(Heq1*10**3)
    vHeq2.append(Heq2*10**3)
    
plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,Pmax1)
#plt.ylim(0,15)
plt.plot(vPmax, vHeq1, 'o', fillstyle='none', linewidth=1.5, label='R=5·10^12 kg·m^-4/s')
plt.plot(vPmax, vHeq2, 'o', linewidth=1.5, label='R=1·10^14 kg·m^-4/s')
plt.xlabel('Pmax rel (mbar)')
plt.ylabel('Heq (mm)')
plt.legend(loc = 1)
plt.show()

