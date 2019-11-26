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
Tin  = 318 #K
Tout = 248 #K
D    = 5.5e-4 #m
massFlow = 2.0 #kg/h

massFlow = massFlow/3600.0 #kg/s
A = math.pi * D * D * 0.25
G = massFlow / A
N = 20

Pin  = PropsSI('P','T',Tin ,'Q',0,fluid)
Pout = PropsSI('P','T',Tout,'Q',0,fluid) #P é cte para qualquer valor de Q(titulo)

deltaP = (Pout - Pin)/N

P1 = Pin
P2 = P1 + deltaP
x_1 = 0

x_2   = getGasFraction(fluid,G,P1,P2,x_1)

rho_1 = getRho(fluid,P1,x_1)
rho_2 = getRho(fluid,P2,x_2)
mu_1  = getMu( fluid,P1,x_1)
mu_2  = getMu( fluid,P2,x_2)
V_1 = G/rho_1
V_2 = G/rho_2

Re_1 = (G*D)/mu_1
Re_2 = (G*D)/mu_2

f_1  = 0.3316/((Re_1)**0.25)
f_2  = 0.3316/((Re_2)**0.25)

f_m = 0.5 * (f_1 + f_2)
V_m = 0.5 * (V_1 + V_2)
deltaL = (deltaP + G*(V_1-V_2))*D/(0.5*f_m*G*V_m)
print(deltaL)
