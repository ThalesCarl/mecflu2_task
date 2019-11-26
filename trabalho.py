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
from CoolProp.CoolProp import PropsSI
from auxiliar import getGasFraction, getRho, getMu

#input data
fluid= 'R600a'
T_in  = 318.15   #K
T_out = 248.15   #K
D     = 5.5e-4 #m
massFlow = 2.0 #kg/h

massFlow = massFlow/3600.0 #kg/s
A = math.pi * D * D * 0.25
G = massFlow / A
N = 70

Pin  = PropsSI('P','T',T_in ,'Q',0,fluid)
Pout = PropsSI('P','T',T_out,'Q',0,fluid) 
deltaP = (Pout - Pin)/N

P_1 = Pin
deltaL = 1.0
while(deltaL > 0 and stopCounter < 1000)
    P_2 = P_1 + deltaP
    x_1 = 0
    
    x_2   = getGasFraction(fluid,G,P_1,P_2,x_1)
    
    rho_1 = getRho(x_1,fluid,P_1)
    rho_2 = getRho(x_2,fluid,P_2)
    mu_1  = getMu(x_1,fluid,P_1)
    mu_2  = getMu(x_2,fluid,P_2)
    V_1 = G/rho_1
    V_2 = G/rho_2
    
    Re_1 = (rho_1*V_1*D)/mu_1
    Re_2 = (rho_2*V_2*D)/mu_2
    
    f_1  = 0.3316/((Re_1)**0.25)
    f_2  = 0.3316/((Re_2)**0.25)
    
    f_m = 0.5 * (f_1 + f_2)
    V_m = 0.5 * (V_1 + V_2)
    
    deltaL = (P_1 - P_2 + G*(V_1-V_2))*D/(0.5*f_m*G*V_m)
print(deltaL)
