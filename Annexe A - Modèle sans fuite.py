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

Q1=10*(10**-9)/60 # m^3/s     Débit 1
Q2=25*(10**-9)/60 # m^3/s     Débit 2
ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique

# Pour les intégrations

temps=10000 #s                Intervalle de temps évalué
pas = 0.1 #s                    Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Définition des EDO's 

F1 = lambda t, H: Q1/A-Patm/(R*A)*(1-ho/Hp)*(1/(1-H/Hp)-1/(1+H/Hp-2*ho/Hp))-2*ro*g/(R*A)*(H-ho)
F2 = lambda t, H: Q2/A-Patm/(R*A)*(1-ho/Hp)*(1/(1-H/Hp)-1/(1+H/Hp-2*ho/Hp))-2*ro*g/(R*A)*(H-ho)

# Résolution des EDO's en utilisant le solveur 
# Méthode Runge-Kutta implicite ordre 5


t_eval = np.arange(0, temps, pas)
sol1 = solve_ivp(F1, [0, temps], [ho], method='Radau', t_eval=t_eval, dense_output=True)
sol2 = solve_ivp(F2, [0, temps], [ho], method='Radau', t_eval=t_eval, dense_output=True)


# Tracer les courbes de H selon t

plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,temps)
#plt.ylim(0.1, 0.2)
plt.plot(sol1.t,  (sol1.y[0])*1000, linewidth=1.5, label='Q1=10 microL/min')
plt.plot(sol2.t, (sol2.y[0])*1000, linestyle = 'dashdot', color = '0', linewidth=1.5, label='Q2=25 microL/min')
plt.xlabel('t (s)')
plt.ylabel('H (mm)')
plt.legend(loc = 4)
plt.show()

#_________________________________________________________________________________________

# Delta P

# Calculer les valeurs de DeltaP dans l'instant final

deltaP1 = Patm*(1-ho/Hp)*(1/(1-sol1.y/Hp)-1/(1+sol1.y/Hp-2*ho/Hp))+2*ro*g*(sol1.y-ho)
deltaP2 = Patm*(1-ho/Hp)*(1/(1-sol2.y/Hp)-1/(1+sol2.y/Hp-2*ho/Hp))+2*ro*g*(sol2.y-ho)

# Tracer les courbes de DeltaP selon t


plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,temps)
#plt.ylim(0,60)
plt.plot(sol1.t, deltaP1[0]*10**-2, label='Q1=10 microL/min')
plt.plot(sol2.t, deltaP2[0]*10**-2, linestyle = 'dashdot', color = '0', linewidth=1.5, label='Q2=25 microL/min')
plt.xlabel('t (s)')
plt.ylabel('Delta P (mbar)')
plt.legend(loc = 4)
plt.show()

#_____________________________________________________


a=2*ro*g
b=-R*Q2-6*ro*g*ho
c=2*R*Q2*ho-2*Patm*(Hp-ho)-2*ro*g*((Hp**2)-2*ho*Hp-2*(ho**2))
d=R*Q2*Hp*(Hp-2*ho)+2*Patm*ho*(Hp-ho)-2*ro*g*Hp*ho*(2*ho-Hp)

f1 = lambda Hi: a*Hi**3 + b*Hi**2 + c*Hi + d

Heq1 = fsolve(f1, Hp)

H = np.linspace(1/1000,11/1000)

plt.figure(figsize = (10,6))
plt.grid()
#plt.xlim(0,temps)
#plt.ylim(0,2300)
plt.plot(H, a*H**3 + b*H**2 + c*H + d, linestyle='solid', color='0', linewidth=1.5)
plt.xlabel('H (m)')
plt.ylabel('dH/dt (m/s)')
plt.show()

print(Heq1*1000)
print(sol2.y[0,n1]*1000)
