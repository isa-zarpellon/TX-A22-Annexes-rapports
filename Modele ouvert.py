import matplotlib.pyplot as plt
import numpy as np
import math


# Données

A=2.01*10**-4 # m^2           Aire de la section transversale
ho=9.95*10**-3 # m            Hauteur intiale de la colomne de liquide dans les puits
Hp=16.5*10**-3 # m            Hauteur des puits
R=10*10**13 # Pa*s/m^3         Résistance hydraulique (considérée constante)

Q1=10*(10**-9)/60 # m^3/s     Débit 1
Q2=25*(10**-9)/60 # m^3/s     Débit 2
ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique
temps=10**7

t = np.linspace(0,10000000)

plt.figure(figsize = (10,6))
plt.grid()
plt.xlim(0,temps)
#plt.ylim(0,0.6)
plt.plot(t, 10**3*(ho+R*Q2/(2*ro*g)-R*Q2/(2*ro*g)*np.exp(-t*2*ro*g/(R*A))), linestyle='solid', color='0', linewidth=1.5)
plt.xlabel('t (s)')
plt.legend(loc = 2)
plt.show()
