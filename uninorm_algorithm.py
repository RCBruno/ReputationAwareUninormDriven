# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 21:15:54 2024

@author: bruno
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import math

def rep_coeff(uv, tv):
    if uv == 0.0:
        return 1.0
    else:
        return 1.0 - uv / tv

def generate_coeff(num_coeff, a, b):
    cont = 0
    coeffs = np.zeros(num_coeff)
    for i in range(1, num_coeff + 1):
        #print("prueba:")
        num_alea = random.randrange(a, b)
        if num_alea % 2 == 1 and num_alea % 3 == 0:
            cont += 1
            coeffs[i-1] = 1.0 - float(cont) / float(i)
        if cont != 0:
            coeffs[i-1] = 1.0 - float(cont) / float(i)
        else:
            coeffs[i-1] = 1.0
    return coeffs

def Ufodor(x, y, neutro):
    #print("RepCoefficient:", x, "RepWeight:", y)
    if 0.0 <= x <= neutro and 0.0 <= y <= neutro:
        return float(2*x*y)
    elif neutro <= x <= 1.0 and neutro <= y <= 1.0:
        return float(2.0*(x+y-x*y-neutro))
    else:
        return float((x+y)/2.0)

def Ucross(x, y):
    if (x == 0.0 and y == 1.0) or (x == 1.0 and y == 0.0):
        return 0.0
    else:
        return (x*y)/(x*y+(1.0-x)*(1.0-y))

def UcrossFodor(x, y):
    if (x == 0.0 and y == 1.0) or (x == 1.0 and y == 0.0):
        return 0.0
    else:
        return (0.5)*(x*y)/(0.5*x*y+0.5*(1.0-x)*(1.0-y))

def U(x, y, neutro):
    print("RepCoefficient:", x, "RepWeight:", y)
    if 0.0 <= x <= neutro and 0.0 <= y <= neutro:
        return float(x*y)
    elif neutro <= x <= 1.0 and neutro <= y <= 1.0:
        return float(x+y-2.0*(x*y))
    else:
        return float((x+y)/2.0)

def rep_weight(neutro, repcoefficient, repweightbefore, round):
    if round == 1.0:
        return neutro
    else:
        return Ufodor(repcoefficient, repweightbefore, neutro)

def rep_weight1(repcoefficient, repweightbefore, round):
    if round == 1.0:
        return 1.0
    else:
        return UcrossFodor(repcoefficient, repweightbefore)

def SuccVR(numUV, j):
    if numUV == 0:
        return 1
    else:
        return 1 - float(numUV/(j+1))

def FuncMu(CoeffRate,N):
    VectorMu = np.zeros(N)
    cte =  1/(math.e-1)
    for i in range(0,N):
        x = CoeffRate[i]
        if x < 0:
            VectorMu[i] = 0
        elif 0 <= x <= 1.0:
            VectorMu[i] =cte*(math.exp(x)-1)
        else:
            VectorMu[i] = 1
    return VectorMu


def FuncNu(CoeffRate,N):
    VectorNu = np.zeros(N)
    for i in range(0,N):
        x = CoeffRate[i]
        if x < 0:
            VectorNu[i] = 1
        elif 0 <= x <= 1.0:
            VectorNu[i] = math.pow(x-1, 2)
        else:
            VectorNu[i] = 0
    return VectorNu

def FuncGrade(Mu, Nu, alpha,N):
    VectorGrade = np.zeros(N)
    for i in range(0,N):
        VectorGrade[i] = math.pow(Mu[i], alpha)/(Mu[i] + Nu[i])
    return VectorGrade


#validator = [1, 0, 1, 1, 1, 1, 1, 0, 1, 1 ]
validator = [1, 0, 0, 1, 1, 1, 1, 1, 1, 1 ]
NumRound = len(validator)
Coeff_Rrate = np.zeros(NumRound)
numUV = 0
x_aux = np.zeros(NumRound)

for j in range(0,NumRound):
    if validator[j]==0:
        numUV = numUV+1
    Coeff_Rrate[j] = SuccVR(numUV, j)

print("List of coefficients of validator:", Coeff_Rrate)

repbefore = np.zeros(NumRound)
repbefore[0] = 1.0
for j in range(0,NumRound):
    repbefore[j] = rep_weight(0.5, Coeff_Rrate[j], repbefore[j-1], j+1)
    #print("RepFor:", repbefore[j])
print("Reputation weight:", repbefore)



Coeff_Mu = FuncMu(Coeff_Rrate, NumRound)
print("Mu:", Coeff_Mu)
Coeff_Nu = FuncNu(Coeff_Rrate, NumRound)
print("Nu:", Coeff_Nu)

alpha = 1.1
alpha12 = 1.3
alpha13 = 1.5
alpha14 = 1.7
alpha15 = 2.0
Coeff_Grade = FuncGrade(Coeff_Mu, Coeff_Nu, alpha, NumRound)
Coeff_Grade12 = FuncGrade(Coeff_Mu, Coeff_Nu, alpha12, NumRound)
Coeff_Grade13 = FuncGrade(Coeff_Mu, Coeff_Nu, alpha13, NumRound)
Coeff_Grade14 = FuncGrade(Coeff_Mu, Coeff_Nu, alpha14, NumRound)
Coeff_Grade15 = FuncGrade(Coeff_Mu, Coeff_Nu, alpha15, NumRound)

#fodor con nu y mu
repFodor = np.zeros(NumRound)
repFodor[0] = 1.0
for j in range(0,NumRound):
    repFodor[j] = rep_weight(0.5, Coeff_Grade[j], repFodor[j-1], j+1)
print("Reputation weight:", repFodor)

repFodor12 = np.zeros(NumRound)
repFodor12[0] = 1.0
for j in range(0,NumRound):
    repFodor12[j] = rep_weight(0.5, Coeff_Grade12[j], repFodor12[j-1], j+1)
print("Reputation weight:", repFodor12)

repFodor13 = np.zeros(NumRound)
repFodor13[0] = 1.0
for j in range(0,NumRound):
    repFodor13[j] = rep_weight(0.5, Coeff_Grade13[j], repFodor13[j-1], j+1)
print("Reputation weight:", repFodor13)

repFodor14 = np.zeros(NumRound)
repFodor14[0] = 1.0
for j in range(0,NumRound):
    repFodor14[j] = rep_weight(0.5, Coeff_Grade14[j], repFodor14[j-1], j+1)
print("Reputation weight:", repFodor14)


repFodor15 = np.zeros(NumRound)
repFodor15[0] = 1.0
for j in range(0,NumRound):
    repFodor15[j] = rep_weight(0.5, Coeff_Grade15[j], repFodor15[j-1], j+1)
print("Reputation weight:", repFodor15)

# Create auxiliary vector to plot
for j in range(0,NumRound):
    x_aux[j]= j+1



plt.plot(x_aux,Coeff_Rrate, '--', label = "$SuccVR(t_1,j)$", color="purple")
#plt.title('Reputation Degree vs. Reputation Coefficient')
#plt.plot(x_aux, repbefore, '-', label = "$Rweight$")
#plt.plot(x_aux, Coeff_Mu, '.-', label = "$Mu$")
#plt.plot(x_aux, Coeff_Nu, ':', label = "$Nu$")
plt.plot(x_aux, Coeff_Grade, ':', label = "$RepD(t_1,j)$", color="brown")
plt.plot(x_aux, repFodor, '.-', label = "$w(t_1,j), alpha=1.1$", color="red")
plt.plot(x_aux, repFodor12, '.-', label = "$w(t_1,j), alpha=1.3$")
plt.plot(x_aux, repFodor13, '.-', label = "$w(t_1,j), alpha=1.5$")
plt.plot(x_aux, repFodor14, '.-', label = "$w(t_1,j), alpha=1.7$")
plt.plot(x_aux, repFodor15, '.-', label = "$w(t_1,j), alpha=2.0$")
plt.xlabel('Round(j)')
#plt.ylabel('Reputation Weight')
plt.grid(True)
even_ticks = [i for i in range(1, NumRound+1) if i % 1 == 0]
plt.xticks(even_ticks)
plt.legend()
plt.savefig("t3.png", dpi=300, bbox_inches='tight', pad_inches=.1)
plt.show()
