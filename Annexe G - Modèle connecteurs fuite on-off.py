import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

plt.style.use('grayscale')

# Données

A=2.01*10**-4 # m^2           Aire de la section transversale
ho=9.95*10**-3 # m            Hauteur initiale de la colonne de liquide dans les puits
Hp=16.5*10**-3 # m            Hauteur des puits
R=10*10**13 # Pa*s/m^3         Résistance hydraulique (considérée constante)

Q2=25*(10**-9)/60 # m^3/s     Débit 2
ro=1000 # kg/m^3              Masse volumique du fluide
g=9.8 # m/s^2                 Accélération de la pesanteur
Patm=101325 # Pa              Pression atmosphérique
 

# Pour les intégrations

temps=10200 #s                Intervalle de temps évalué
pas = 1 #s                    Pas
n=int(temps/pas) #            Nombre d'instants considérés
n1=int(n-1) #                 n-1

# Calcul de Heq

a=2*ro*g
b=-R*Q2-6*ro*g*ho
c=2*R*Q2*ho-2*Patm*(Hp-ho)-2*ro*g*(Hp**2-2*ho*Hp-2*ho**2)
d=R*Q2*Hp*(Hp-2*ho)+2*Patm*ho*(Hp-ho)+2*ro*g*ho*Hp*(Hp-2*ho)

f = lambda Hi: a*Hi**3 + b*Hi**2 + c*Hi + d

Heq = fsolve(f, Hp)

# Def de la fonction "fuite"

def fuite (q, Pseuil) :

    C=2*ho
    
    F1 = lambda t, H: Q2/A-Patm/(R*A)*(Hp-ho)*(1/(Hp-H)-1/(Hp+H-C))-ro*g/(R*A)*(2*H-C)
    F2 = lambda t, H: Q2/A-q/A-Patm/(R*A)*(Hp-ho)*(1/(Hp-H)-1/(Hp+H+q*(t-t0)/A-C))-ro*g/(R*A)*(2*H-C+q*(t-t0)/A)
    
    
    t_evala = np.arange(0, temps, pas)
    sol = solve_ivp(F1, [0, temps], [ho], method='Radau', t_eval=t_evala)
    
    model_der = "f1"
    Pinf = 0
    j_der = 0
    i = 0
    sol = solve_ivp(F1, [0, temps], [ho], method='Radau', t_eval=t_evala)
    h_values = [ho/Heq]
    Q_prime = [0]
    H_total = [0]
    
    for j in range (0, n1, 1) :
        if C >= 0 :
            if i == 0 : 
                Pinf=ro*g*sol.y[0,j]+Patm*(Hp-ho)/(Hp-sol.y[0,j-j_der]) 
                i += 1
                
            elif Pinf < Pseuil :
               
                if model_der == "f1" : 
                    Pinf=ro*g*sol.y[0,j-j_der]+Patm*(Hp-ho)/(Hp-sol.y[0,j-j_der])
                    h_values.append(sol.y[0,j-j_der]/Heq)
                    Qp=Patm/R*(Hp-ho)*(1/(Hp-sol.y[0,j-j_der])-1/(Hp+sol.y[0,j-j_der]-C))+ro*g/R*(2*sol.y[0,j-j_der]-C)
                    Q_prime.append(Qp/Q2)    
                    H_total.append(C/(2*ho))
                else: 
                    # print("F1 : " + str(j*pas))
                    model_der = "f1"
                    C=C-q/A*(j-j_der)*pas
                    t_evala = np.arange(j*pas, temps, pas)
                    sol = solve_ivp(F1, [j*pas, temps+1], [sol.y[0,j-j_der]], method='Radau', t_eval=t_evala)
                    j_der = j - 1
                    Pinf=ro*g*sol.y[0,j-j_der]+Patm*(Hp-ho)/(Hp-sol.y[0,j-j_der]) 
                    h_values.append(sol.y[0,j-j_der]/Heq)
                    Qp=Patm/R*(Hp-ho)*(1/(Hp-sol.y[0,j-j_der])-1/(Hp+sol.y[0,j-j_der]-C))+ro*g/R*(2*sol.y[0,j-j_der]-C)
                    Q_prime.append(Qp/Q2) 
                    H_total.append(C/(2*ho))
            else :
                if model_der == "f2" : 
                    Pinf=ro*g*sol.y[0,j-j_der]+Patm*(Hp-ho)/(Hp-sol.y[0,j-j_der])
                    h_values.append(sol.y[0,j-j_der]/Heq)
                    Qp=Patm/R*(Hp-ho)*(1/(Hp-sol.y[0,j-j_der])-1/(Hp+sol.y[0,j-j_der]+q*(j-j_der)*pas/A-C))+ro*g/R*(2*sol.y[0,j-j_der]-C+q*(j-j_der)*pas/A)
                    Q_prime.append(Qp/Q2) 
                    H_total.append((C-q*(j-j_der)*pas/A)/(2*ho))
                else :
                    # print("F2 : " + str(j*pas))
                    model_der = "f2"
                    t0=j*pas
                    t_evala = np.arange(j*pas, temps, pas)
                    #ho = sol.y[0,j-j_der]
                    sol = solve_ivp(F2, [j*pas, temps+1], [sol.y[0,j-j_der]], method='Radau', t_eval=t_evala)
                    j_der = j - 1
                    Pinf=ro*g*sol.y[0,j-j_der]+Patm*(Hp-ho)/(Hp-sol.y[0,j-j_der]) 
                    h_values.append(sol.y[0,j-j_der]/Heq)
                    Qp=Patm/R*(Hp-ho)*(1/(Hp-sol.y[0,j-j_der])-1/(Hp+sol.y[0,j-j_der]+q*(j-j_der)*pas/A-C))+ro*g/R*(2*sol.y[0,j-j_der]-C+q*(j-j_der)*pas/A)
                    Q_prime.append(Qp/Q2) 
                    H_total.append((C-q*(j-j_der)*pas/A)/(2*ho))
        else :
            print("Vol total = 0" + str(j*pas))
            break
    
    return (h_values, Q_prime, H_total)
    
# Calcul des Hauteurteurs et des débits 

sol11, Qp11, Ht11 = fuite(0.01*Q2, 100*100+Patm)
sol21, Qp21, Ht21 = fuite(0.05*Q2, 100*100+Patm)
sol31, Qp31, Ht31 = fuite(0.10*Q2, 100*100+Patm)
sol41, Qp41, Ht41 = fuite(0.50*Q2, 100*100+Patm)

sol12, Qp12, Ht12 = fuite(0.01*Q2, 200*100+Patm)
sol22, Qp22, Ht22 = fuite(0.05*Q2, 200*100+Patm)
sol32, Qp32, Ht32 = fuite(0.10*Q2, 200*100+Patm)
sol42, Qp42, Ht42 = fuite(0.50*Q2, 200*100+Patm)



# Courbes 


t_v= np.arange(0, temps-pas, pas)

def plot (f1, f2, f3, f4, y, vP, ymin, ymax) :

    plt.figure(figsize = (10,6))
    plt.grid()
    plt.xlim(0,temps)
    plt.ylim(ymin, ymax)
    plt.plot(t_v, f1, color='0', linestyle='dashdot', label='1%')
    plt.plot(t_v, f2, color='0', linestyle='solid', label='5%')
    plt.plot(t_v, f3, color='0', linestyle='dotted', label='10%')
    plt.plot(t_v, f4, color='0', linestyle=(0, (5, 5)), label='50%')
    plt.legend(loc = 4, title='Pression seuil égale à {} mbar - q/Q :'.format(vP))
    plt.ylabel('{}'.format(y))
    plt.xlabel("temps(s)")
    plt.show()
    
#plot (f1, f2, f3, f4, titre axe y, Pseuil, ymin axe, ymax axe) :
   
plot(sol11, sol21, sol31, sol41, 'H/Heq', '100 mbar', 0.85, 1.001)
plot(Qp11, Qp21, Qp31, Qp41, "Q'/Q", '100 mbar', 0.8, 1)
plot(Ht11, Ht21, Ht31, Ht41, 'Htotal/Ht0', '100 mbar', 0.85, 1.001)

plot(sol12, sol22, sol32, sol42, 'H/Heq', '200 mbar', 0.85, 1.001)
plot(Qp12, Qp22, Qp32, Qp42, "Q'/Q", '200 mbar', 0.8, 1)
plot(Ht12, Ht22, Ht32, Ht42, 'Htotal/Ht0', '200 mbar', 0.85, 1.001)
