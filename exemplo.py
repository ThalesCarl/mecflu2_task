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
from auxiliar import *

#input data
fluid= 'R22'
T_in  = 40 + 273.15  #K
T_out = 5  + 273.15 #K
D    = 1.63e-3 #m
massFlow = 0.010 #kg/h

A = math.pi * D * D * 0.25
G = massFlow / A
N = 350

P_in  = PropsSI('P','T',T_in ,'Q',0,fluid)
P_out = PropsSI('P','T',T_out,'Q',0,fluid) 

deltaP = (P_out - P_in)/N
deltaT = (T_out - T_in)/N

T_1 = T_in
x_1 = 0
L_final = 0
for i in range(N):
    T_2 = T_1 + deltaT
    P_1 = PropsSI('P','T',T_1 ,'Q',0,fluid)
    P_2 = PropsSI('P','T',T_2 ,'Q',0,fluid)
    
    x_2   = getGasFractionUsingPressure(fluid,G,P_1,P_2,x_1)    
    
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
    
    
    
    x_1 = x_2
    T_1 = T_2
    L_final += deltaL
print(L_final)
    
