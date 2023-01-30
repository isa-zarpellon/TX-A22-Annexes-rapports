import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

plt.style.use('grayscale')

# Données

A=2.01*10**-4 # m^2           Aire de la section transversale
ho=9.95*10**-3 # m            Hauteur initiale de la colonne de liquide dans les puits
Hp=16.5*10**-3 # m            Hauteur des puits
R=10*10**13 # Pa*s/m^3        Résistance hydraulique (considérée constante)

Q1=10*(10**-9)/60 # m^3/s     Débit 1
Q2=10*(10**-9)/60 # m^3/s     Débit 2
ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique
Pmax=50*100+Patm # Pa         Pression maximale dans le puits d'entrée
Pmin=-100*100+Patm # Pa       Dépression maximal dans le puits de sortie

# Pour les intégrations

temps=10000 #s                Intervalle de temps évalué
pas = 0.01 #s                 Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Définition des EDO's 

F1 = lambda t, H: Q1/A-Patm/(R*A)*(1-ho/Hp)*(1/(1-H/Hp)-1/(1+H/Hp-2*ho/Hp))-2*ro*g/(R*A)*(H-ho)
F2 = lambda t, H: Q2/A-Patm/(R*A)*(1-ho/Hp)*(1/(1-H/Hp)-1/(1+H/Hp-2*ho/Hp))-2*ro*g/(R*A)*(H-ho)

# Résolution des EDO's en utilisant le solver 
# Méthode Runge-Kutta implicite ordre 5


t_eval = np.arange(0, temps, pas)
sol1 = solve_ivp(F1, [0, temps], [ho], method='Radau', t_eval=t_eval, dense_output=True)
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



#__________________________________________________________________________________________

# Modèle avec fuite dans les deux puits - résolu pour Q2

# Détermination de l'instant où Ps <= Pmax

k=0

for l in range (0, n1, 1) :
    
    if Patm*(Hp-ho)/(Hp-2*ho+(sol2b.y[0,l])) > Pmin :
        k=k+1
        
    else :
        break

# EDO qui définit le troisième régime

F23 = lambda t, H: Q2/A-(Pmax-Pmin)/(R*A)-2*ro*g/(R*A)*(H-ho)

# Les nouveaux intervalles de temps

t_eval23a = np.arange(0, j*pas, pas)
t_eval23b = np.arange(j*pas, j*pas+k*pas, pas)
t_eval23c = np.arange(j*pas+k*pas, temps, pas)

print(j*pas+k*pas)
print(sol2b.y[0,k])

# Recalcul de la première solution avec le nouvel intervalle 
# Cacul de la deuxième solution

sol23a = solve_ivp(F2, [0, j*pas], [ho], method='Radau', t_eval=t_eval23a)
sol23b = solve_ivp(F22, [j*pas, j*pas+k*pas+1], [sol2.y[0,j]], method='Radau', t_eval=t_eval23b)
sol23c = solve_ivp(F23, [j*pas+k*pas, temps+1], [sol2b.y[0,k]], method='Radau', t_eval=t_eval23c)

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
plt.plot(sol23a.t, (sol23a.y[0])*1000, linestyle='dotted', color='0', linewidth=3, label="Modèle avec fuite dans les deux puits")
plt.plot(sol23b.t, (sol23b.y[0])*1000, linestyle='dotted', color='0', linewidth=3)
plt.plot(sol23c.t, (sol23c.y[0])*1000, linestyle='dotted', color='0', linewidth=3)
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
plt.plot(sol23a.t, (Patm*(Hp-ho)/(Hp+sol23a.y[0]-2*ho)-Patm)*10**-2, linestyle='dotted', color='0', linewidth=3, label="Modèle avec fuite dans les deux puits")
plt.plot(sol23b.t, (Patm*(Hp-ho)/(Hp+sol23b.y[0]-2*ho)-Patm)*10**-2, linestyle='dotted', color='0', linewidth=3)
plt.plot(sol23c.t, (Pmin-Patm)*(10**-2)*np.ones(np.shape(sol23c.t)), linestyle='dotted', color='0', linewidth=3)
plt.xlabel('t (s)')
plt.ylabel('Ps rel (mbar)')
plt.legend(loc = 1)
plt.show()
