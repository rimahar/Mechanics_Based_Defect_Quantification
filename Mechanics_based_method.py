# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:45:54 2023

@author: rsa87
"""
import numpy as np
from scipy.interpolate import LinearNDInterpolator
outputFile = 'Data_All.csv'
results=np.genfromtxt(outputFile, delimiter=",")
E_11=results[:,3] #E11
E_33=results[:,4] #E33
nu_13=results[:,5] #nu13
voids=results[:,0] #Total Volume of Voids
theta=np.multiply(results[:,1],2) #Alignment
pores_voids=results[:,2]  #Pores
Voids=LinearNDInterpolator((E_11, E_33, nu_13), voids)
Theta_t=LinearNDInterpolator((E_11, E_33, nu_13), theta)
Pores_Voids_Ratio=LinearNDInterpolator((E_11,E_33, nu_13), pores_voids)


