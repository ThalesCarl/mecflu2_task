#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Aluno: Thales Carl Lavoratti
Matricula: 15100656
Data: 27/11/2019

Trabalho computacional da disciplina EMC5419 – MECÂNICA DOS FLUIDOS II
"""
#import numpy as 
import math
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from auxiliar import getGasFraction, getRho, getMu, getS, getH, SolutionContext

#input data
fluid= 'R600a'
T_in  = 318.15   #K
T_out = 248.15   #K
D     = 5.5e-4 #m
massFlow = 2.0 #kg/h

N = 45
massFlow = massFlow/3600.0 #kg/s
A = math.pi * D * D * 0.25
G = massFlow / A

P_in  = PropsSI('P','T',T_in ,'Q',0,fluid)
P_out = PropsSI('P','T',T_out,'Q',0,fluid) 
deltaP = (P_out - P_in)/N


P_1 = P_in
deltaL = 1.0
L = 0
x_1 = 0
maxIterations = 10000
stopCounter = 0
save = SolutionContext()
while(deltaL > 0):
    P_2 = P_1 + deltaP
        
    x_2   = getGasFraction(fluid,G,P_1,P_2,x_1)
    rho_1 = getRho(x_1,fluid,P_1)
    rho_2 = getRho(x_2,fluid,P_2)
    mu_1  = getMu(x_1,fluid,P_1)
    mu_2  = getMu(x_2,fluid,P_2)
    s_1   = getS(x_1,fluid,P_1)
    h_1   = getH(x_1,fluid,P_1)
    save.s.append(s_1)
    save.h.append(h_1)
    
    V_1 = G/rho_1
    V_2 = G/rho_2
    
    Re_1 = (rho_1*V_1*D)/mu_1
    Re_2 = (rho_2*V_2*D)/mu_2
    
    f_1  = 0.3316/((Re_1)**0.25)
    f_2  = 0.3316/((Re_2)**0.25)
    
    f_m = 0.5 * (f_1 + f_2)
    V_m = 0.5 * (V_1 + V_2)
    
    deltaL = (P_1 - P_2 + G*(V_1-V_2))*D/(0.5*f_m*G*V_m)
    L += deltaL
    
    x_1 = x_2
    P_1 = P_2
    stopCounter += 1
    if(stopCounter > maxIterations):
        print('Código não convergiu')
        deltaL = float('nan')
        
  
plt.grid(color='k', linestyle=':', linewidth=0.5)
plt.plot(save.s,save.h,'ko-')
plt.xlabel('s [J/kgK]')
plt.ylabel(r'h [J/kg]')
plt.savefig('mollier.png')
print(L)
print(s_2)
