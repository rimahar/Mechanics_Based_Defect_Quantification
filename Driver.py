# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 23:39:57 2024

@author: RIMAH
"""
from MechanicsModel import MechanicsModel
#Material Propeties
E0 = 70000
nu0 = 0.33
#Inlcusion Properties
E1 = 0
nu1 = 0
#Voids Volume Fraction
VoidsVolumeFraction_List = [0.01,0.05]
#Pores Fraction of total Voids popullation
PoresFraction_List = [0.1]
#ODF
Theta_List = [81]
##Aspect Ratio of Incusion
alpha_List = [100]
##Run the model
data_processor = MechanicsModel(E0=E0, nu0=nu0, E1=E1, nu1=nu1, VoidsVolumeFraction_List=VoidsVolumeFraction_List, 
  PoresFraction_List=PoresFraction_List, Theta_all=Theta_List, alpha_List=alpha_List, path='D:/SRNL/Analytical Model/Testing/',
  NUMBER_OF_INCLUSIONS=2, OUTPUT_FILE='Data_2.csv')
data_processor.run_Model()