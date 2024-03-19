# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:50:05 2024

@author: RIMAH
"""

import os
import csv

import numpy as np
import numpy as np
from numpy import arccosh
from numpy import arccos
from numpy import linspace
from numpy import pi
from numpy import logspace
from numpy import gradient
import scipy.optimize
import scipy.integrate as integrate
import math
class MechanicsModel:
    def __init__(self, E0, nu0, E1, nu1, VoidsVolumeFraction_List, PoresFraction_List, Theta_all, alpha_List, path='D:/SRNL/Analytical Model/Testing/', NUMBER_OF_INCLUSIONS=2, OUTPUT_FILE='Data_1.csv'):
        self.E0 = E0
        self.nu0 = nu0
        self.E1 = E1
        self.nu1 = nu1
        self.path = path
        self.NUMBER_OF_INCLUSIONS = NUMBER_OF_INCLUSIONS
        self.OUTPUT_FILE = OUTPUT_FILE
        self.VoidsVolumeFraction_List = VoidsVolumeFraction_List
        self.PoresFraction_List = PoresFraction_List
        self.Theta_all = Theta_all
        self.alpha_List = alpha_List
        self.headers = ['VoidsVolumeFraction', 'PoresFraction', 'Theta_t', 'Alpha', 'E_11', 'E_33','nu_13', 'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8', 'C_9']

    def run_Model(self):
        try:
            os.chdir(self.path)

            # Assuming KK is a tuple or list returned by Dr.Run
            KK = self.Run(self.E0, self.nu0, self.E1, self.nu1, self.VoidsVolumeFraction_List, self.PoresFraction_List, self.Theta_all, self.alpha_List, self.NUMBER_OF_INCLUSIONS + 1)
           # print(KK)
            # Transpose the data
         #   outputData = list(map(list, zip(self.VoidsVolumeFraction_List, self.PoresFraction_List, self.Qall, self.alpha_List, *KK)))
           # outputData = np.transpose([self.VoidsVolumeFraction_List, self.PoresFraction_List, self.Qall, self.alpha_List])
            outputData=np.transpose([KK[0],KK[1],KK[2],KK[3],KK[4],KK[5],KK[6],KK[7],KK[8],KK[9],KK[10],KK[11],KK[12],KK[13],KK[14],KK[15]])
            print(outputData)
           # print(outputData2)

            try:
                results = np.genfromtxt(self.OUTPUT_FILE, delimiter=",")
                with open(self.OUTPUT_FILE, mode='a', newline='') as csvfile:
                    output = csv.writer(csvfile, delimiter=',', quotechar='|')

                    # Write headers
                  #  output.writerow(self.headers)

                    # Write data
                    for row in outputData:
                        output.writerow(row)
            except:
                with open(self.OUTPUT_FILE, mode='w', newline='') as csvfile:
                    output = csv.writer(csvfile, delimiter=',', quotechar='|')

                    # Write headers
                    output.writerow(self.headers)

                    # Write data
                    for row in outputData:
                        output.writerow(row)
        
            # Writing to CSV file with headers
            

        except Exception as e:
            print(f"An error occurred: {e}")
    
    
            
    def Run(self,E0,nu0,E1,nu1,VoidsVolumeFraction_List,PoresFraction_List,Theta_all,alpha_List,numberOfinclusions):
        jj=0 
        E111=[]
        E222=[]
        C00=[]    
        C11=[]    
        C22=[]    
        C33=[]    
        C44=[]    
        C55=[]    
        C66=[]    
        C77=[]    
        C88=[]
        nu21_List=[]
        nu31_List=[]
        nu32_List=[]
        nu12_List=[]
        nu13_List=[]
        nu23_List=[]
        VoidsVolumeFractions=[]
        PoresFractions=[]
        Qs=[]
        CrackAspectRatios=[]
        for VoidsVolumeFraction in VoidsVolumeFraction_List:
            for Theta_t in Theta_all:
                    for PoresFraction in PoresFraction_List:
                        for CrackAspectRatio in alpha_List:
                            print('-----------------------------------')
                            print('Voids Volume Fraction: '+str(VoidsVolumeFraction))
                            print('Pores/Voids: '+str(PoresFraction))
                            print('Q: '+str(Theta_t))
                            print('Aspect Ratio: '+str(CrackAspectRatio))
                            print('-----------------------------------')
                            Q=self.getQ(Theta_t)
                            C=self.ODF(E0,nu0,E1,nu1,VoidsVolumeFraction,PoresFraction,Q,numberOfinclusions,CrackAspectRatio)    
                            
                            f1 = -C[1] * C[4]**2. + 2. * C[3] * C[4] * C[5] - C[0] * C[5]**2. + C[2] * (C[0] * C[1] - C[3]**2.)
                            f2 = (C[4] * C[5] - C[3] * C[2]) / f1
                            f3 = (C[3] * C[5] - C[4] * C[1]) / f1
                            f4 = (C[3] * C[4] - C[0] * C[5]) / f1
                
                            # Set values
                            S = np.empty(9)
                            S[0] = (C[1] * C[2] - C[5]**2.) / f1
                            S[1] = (C[0] * C[2] - C[4]**2.) / f1
                            S[2] = (C[0] * C[1] - C[3]**2.) / f1
                            S[3] = f2
                            S[4] = f3
                            S[5] = f4
                            S[6] = 1. / C[6]
                            S[7] = 1. / C[7]
                            S[8] = 1. / C[8]
                            
                            E = [1./S[i] for i in range(3)]
                            G = [1./S[i] for i in range(6,9)]
                            nu21 = -E[1]*S[3]
                            nu31 = -E[2]*S[4]
                            nu32 = -E[2]*S[5]
                            nu12 = -E[0]*S[3]
                            nu13 = -E[0]*S[4]
                            nu23 = -E[1]*S[5]
                            E11 = 1 / S[2]
                            E22 = 1 / S[0]
                            
                            # Append results to lists
                            C00.append(C[0])
                            C11.append(C[1])
                            C22.append(C[2])
                            C33.append(C[3])
                            C44.append(C[4])
                            C55.append(C[5])
                            C66.append(C[6])
                            C77.append(C[7])
                            C88.append(C[8])
                            
                            nu21_List.append(nu21)
                            nu31_List.append(nu31)
                            nu32_List.append(nu32)
                            nu12_List.append(nu12)
                            nu13_List.append(nu13)
                            nu23_List.append(nu23)
                            
                            
                            E111.append(float(E11))
                            E222.append(float(E22))  
                            VoidsVolumeFractions.append(VoidsVolumeFraction)
                            PoresFractions.append(PoresFraction)
                            Qs.append(Theta_t)
                            CrackAspectRatios.append(CrackAspectRatio)
                            
                            print(E111)
        return  VoidsVolumeFractions, Qs,PoresFractions,CrackAspectRatios,E111,E222,nu13_List,C00,C11,C22,C33,C44,C55,C66,C77,C88       

    def ODF(self,E0,nu0,E1,nu1,c_volume,CPores,Q,numberOfinclusions,alpha1):
        #### Begin User inputs ####
        P = .5

         ##### End User inputs #####
        alpha1=1/alpha1 
        # #Defining aij
        a11=lambda beta, phi, theta:  np.cos(theta)*np.cos(phi)*np.cos(beta)-np.sin(phi)*np.sin(beta)
        a12=lambda beta, phi, theta:  -np.cos(theta)*np.cos(phi)*np.sin(beta)-np.sin(phi)*np.cos(beta)
        a13=lambda beta, phi, theta:  np.sin(theta)*np.cos(phi)
        a21=lambda beta, phi, theta:  np.cos(theta)*np.sin(phi)*np.cos(beta)+np.cos(phi)*np.sin(beta)
        a22=lambda beta, phi, theta:  -np.cos(theta)*np.sin(phi)*np.sin(beta)+np.cos(phi)*np.cos(beta)
        a23=lambda beta, phi, theta:  np.sin(theta)*np.sin(phi)
        a31=lambda beta, phi, theta:  -np.sin(theta)*np.cos(beta)
        a32=lambda beta, phi, theta:  np.sin(theta)*np.sin(beta)
        a33=lambda beta, phi, theta:  np.cos(theta)
        a=[[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]   
        
        #Defining the denomenator in g:    
        g_den=lambda thetaI: np.sin(thetaI/2)**(2*P-1)*np.cos(thetaI/2)**(2*Q-1)#*2
        gg_den=integrate.quad(g_den,0,np.pi/2)  # Limit here needs to be consistent with that in the triple integral
        g=lambda theta:((((np.sin(theta/2))**(2*P-1))*((np.cos(theta/2))**(2*Q-1)))/(gg_den[0]))#/2#3.14/90#(8.*np.pi)
        thetaa=np.arange(0,np.pi/2,0.01)
        ggg=(((np.sin(thetaa/2))**(2*P-1))*((np.cos(thetaa/2))**(2*Q-1)))/(gg_den[0])#*3.14/180#(8.*np.pi)
        thetaa=np.arange(0,np.pi/2,0.01)
        #Defining L_ijkl2
        L_ijkl2=np.zeros(81)
        L_ijkl2=L_ijkl2.reshape(3,3,3,3)
        L_ijkl2_phase2=np.zeros(81)
        L_ijkl2_phase2=L_ijkl2.reshape(3,3,3,3)
        L_ijkl3=np.zeros(81)
        L_ijkl3=L_ijkl3.reshape(3,3,3,3)
        #Defining L_pqrs
        L_pqrs=np.zeros(81)
        L_pqrs=L_pqrs.reshape(3,3,3,3)
        L_matrix = np.zeros((3,3,3,3))
        L_pqrs_phase2=np.zeros(81)
        L_pqrs_phase2=L_pqrs_phase2.reshape(3,3,3,3)
        Cmatrix = self.spheroid(E0, nu0, E1, nu1, alpha=100., c=0.0, type='prolate')
        ## Enter L
        # Unique entries in L
        L_matrix[2][2][2][2] = Cmatrix[1]
        L_matrix[2][2][0][0] = Cmatrix[3]
        L_matrix[0][0][0][0] = Cmatrix[0]
        L_matrix[0][0][1][1] = Cmatrix[2]
        L_matrix[1][1][2][2] = Cmatrix[3]
        L_matrix[1][1][1][1] = Cmatrix[0]
        L_matrix[0][1][0][1] = Cmatrix[5]#/2
        L_matrix[2][1][2][1] = Cmatrix[4]#/2
        L_matrix[2][0][2][0] = Cmatrix[4]#/2
        # Symmetry of off diagonal normal terms
        L_matrix[0][0][2][2]=L_matrix[2][2][0][0]
        L_matrix[2][2][1][1]=L_matrix[1][1][2][2]
        L_matrix[1][1][0][0]=L_matrix[0][0][1][1]
        # Symmetry of shear terms
        L_matrix[0][2][0][2]=L_matrix[2][0][2][0]
        L_matrix[0][2][2][0]=L_matrix[2][0][2][0]
        L_matrix[2][0][0][2]=L_matrix[2][0][2][0]
        L_matrix[1][0][1][0]=L_matrix[0][1][0][1]
        L_matrix[1][0][0][1]=L_matrix[0][1][0][1]
        L_matrix[0][1][1][0]=L_matrix[0][1][0][1]
        L_matrix[1][2][1][2]=L_matrix[2][1][2][1]
        L_matrix[1][2][1][2]=L_matrix[2][1][2][1]
        L_matrix[1][2][2][1]=L_matrix[2][1][2][1]
        L_matrix[2][1][1][2]=L_matrix[2][1][2][1]
        Ccracks=1-CPores
        C=self.model(alpha1,CPores*c_volume,Ccracks*c_volume,numberOfinclusions,E0,nu0,E1,nu1)
        ## Enter L
        # Unique entries in L
        L_pqrs[2][2][2][2] = C[1]
        L_pqrs[2][2][0][0] = C[3]
        L_pqrs[0][0][0][0] = C[0]
        L_pqrs[0][0][1][1] = C[2]
        L_pqrs[1][1][2][2] = C[3]
        L_pqrs[1][1][1][1] = C[0]
        L_pqrs[0][1][0][1] = C[5]
        L_pqrs[2][1][2][1] = C[4]
        L_pqrs[2][0][2][0] = C[4]
        # Symmetry of off diagonal normal terms
        L_pqrs[0][0][2][2]=L_pqrs[2][2][0][0]
        L_pqrs[2][2][1][1]=L_pqrs[1][1][2][2]
        L_pqrs[1][1][0][0]=L_pqrs[0][0][1][1]
        # Symmetry of shear terms
        L_pqrs[0][2][0][2]=L_pqrs[2][0][2][0]
        L_pqrs[0][2][2][0]=L_pqrs[2][0][2][0]
        L_pqrs[2][0][0][2]=L_pqrs[2][0][2][0]
        L_pqrs[1][0][1][0]=L_pqrs[0][1][0][1]
        L_pqrs[1][0][0][1]=L_pqrs[0][1][0][1]
        L_pqrs[0][1][1][0]=L_pqrs[0][1][0][1]
        L_pqrs[1][2][1][2]=L_pqrs[2][1][2][1]
        L_pqrs[1][2][2][1]=L_pqrs[2][1][2][1]
        L_pqrs[2][1][1][2]=L_pqrs[2][1][2][1]     
        LPQRS=[[L_pqrs[0][0][0][0],L_pqrs[0][0][1][1],L_pqrs[0][0][2][2],0,0,0],[L_pqrs[1][1][0][0],L_pqrs[1][1][1][1],L_pqrs[1][1][2][2],0,0,0],[L_pqrs[2][2][0][0],L_pqrs[2][2][1][1],L_pqrs[2][2][2][2],0,0,0],[0,0,0,L_pqrs[1][2][1][2],0,0],[0,0,0,0,L_pqrs[0][2][0][2],0],[0,0,0,0,0,L_pqrs[0][1][0][1]]]
        L_pqrs=L_pqrs-L_matrix
        counter=0
        g2=lambda theta:((((np.sin(theta/2))**(2*P-1))*((np.cos(theta/2))**(2*Q-1)))/(gg_den[0]))#*2#/(3.14/90)
        C=lambda beta,phi, theta: (g2(theta)*np.sin(theta))
        factor=integrate.tplquad(C, 0, np.pi/2, lambda theta: 0, lambda theta: 2*np.pi, lambda theta, phi: 0, lambda theta, phi: 2*np.pi)[0]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        
                        counter = 0
                      
                        if (i==0 and j==0 and k==0 and l==0) or (i==0 and j==0 and k==1 and l==1) or (i==0 and j==0 and k==2 and l==2) or (i==2 and j==2 and k==2 and l==2) or (i==1 and j==2 and k==1 and l==2) or (i==0 and j==1 and k==0 and l==1):
                            for p in range(3):
                                for q in range(3):
                                    for r in range(3):
                                        for s in range(3):
                                            
                                            C=lambda beta,phi, theta:(1*g(theta))*np.sin(theta)*a[i][p](beta,phi,theta)*a[j][q](beta,phi,theta)*a[k][r](beta,phi,theta)*a[l][s](beta,phi,theta)*L_pqrs[p][q][r][s]/factor
                                            C2=integrate.tplquad(C, 0, np.pi/2, lambda theta: 0, lambda theta: 2*np.pi, lambda theta, phi: 0, lambda theta, phi: 2*np.pi)
                    
                                            L_ijkl2[i][j][k][l]=C2[0]+L_ijkl2[i][j][k][l]
                                          
                                            L_ijkl3[i][j][k][l]=1*L_matrix[i][j][k][l]+1*(L_ijkl2[i][j][k][l])
                                         
                                            counter += 1
                                            
                          #  print(str(i)+str(j)+str(k)+str(l))

                          #  print(L_ijkl3[i][j][k][l])
                            
                                            

        
        return([L_ijkl3[0][0][0][0],L_ijkl3[0][0][0][0],L_ijkl3[2][2][2][2],L_ijkl3[0][0][1][1],L_ijkl3[0][0][2][2],L_ijkl3[0][0][2][2],L_ijkl3[1][2][1][2],L_ijkl3[1][2][1][2],L_ijkl3[0][1][0][1]])



    def model(self,alpha,den1,den2,mm,E0,nu0,E1,nu1):
        mu0=nu0
        mu1=nu1
        v=mu0
        o=0
        if o==0: 
            g=(alpha/(1-alpha**2)**(3/2))*(np.cosh(alpha)-alpha*(1-alpha**2)**0.5)
        if o==1:
             g=(alpha/(1-alpha**2)**(3/2))*(math.acos(alpha)+alpha*(1-alpha**2)**0.5)
        
        lambda0 = E0*nu0/((1.+nu0)*(1.-2.*nu0))
        mu0 = E0/(2.*(1.+nu0))
        lambda1 = E1*nu1/((1.+nu1)*(1.-2.*nu1))
        mu1 = E1/(2.*(1.+nu1))
        
        
        
        g = alpha/(1. - alpha**2.)**1.5*(arccos(alpha) - alpha*(1. - alpha**2.)**.5)    
        D1 = 1. + 2.*(mu1 - mu0)/(lambda1 - lambda0)     
        D2 = (lambda0 + 2.*mu0)/(lambda1 - lambda0)     
        D3 = lambda0/(lambda1 - lambda0)     
             
        S1111 = 1./(2.*(1. - nu0))*(1. - 2.*nu0 + (3.*alpha**2. - 1.)/(alpha**2. - 1.) - (1. - 2.*nu0 + 3.*alpha**2./(alpha**2. - 1.))*g)     
             
        S2222 = 3./(8.*(1. - nu0))*alpha**2./(alpha**2. - 1.) + 1./(4.*(1. - nu0))*(1. - 2.*nu0 - 9./(4.*(alpha**2. - 1.)))*g     
        S3333 = S2222     
             
        S2233 = 1./(4.*(1. - nu0))*(alpha**2./(2.*(alpha**2. - 1.)) - (1. - 2.*nu0 + 3./(4.*(alpha**2 - 1.)))*g)     
        S3322 = S2233     
             
        S2211 = -1./(2.*(1. - nu0))*alpha**2./(alpha**2. - 1.) + 1./(4*(1. - nu0))*(3*alpha**2./(alpha**2. - 1.) - (1. - 2.*nu0))*g     
        S3311 = S2211     
             
        S1122 = -1./(2.*(1. - nu0))*(1. - 2.*nu0 + 1./(alpha**2. - 1.)) + 1./(2.*(1. - nu0))*(1. - 2.*nu0 + 3./(2.*(alpha**2. - 1.)))*g     
        S1133 = S1122     
             
        S2323 = 1./(4.*(1.-nu0))*(alpha**2./(2.*(alpha**2. - 1.)) + (1. - 2.*nu0 - 3./(4.*(alpha**2. - 1.)))*g)     
        S3232 = S2323     
             
        S1212 = 1./(4.*(1.-nu0))*(1. - 2.*nu0 - (alpha**2. + 1.)/(alpha**2. - 1.) - .5*(1. - 2.*nu0 - (3.*(alpha**2. + 1.))/(alpha**2. - 1.))*g)     
        S1313 = S1212     
           
        k0=E0/(3*(1-2*v)) 
        S_1111=[1,(7-5*v)/(15*(1-v)),S1111]
     
        S_2222=[1,(7-5*v)/(15*(1-v)),S2222]

        S_3333=[1, (7-5*v)/(15*(1-v)),S3333]#

        S_1122=[1,(5*v-1)/(15*(1-v)),S1122]
        S_1133=[1,0,S1133]#v/(1-v)]
        S_2233=[1,(5*v-1)/(15*(1-v)),S2233]
        S_3322=[1,0, S3322]
        S_2211=[1,0, S2211]
        S_3311=[1, (5*v-1)/(15*(1-v)), S3311]
     
        S_1212=[1,(4-5*v)/(15*(1-v)),S1212]
        S_1313=[ 1,0,S1313]
        S_2323=[1,(4-5*v)/(15*(1-v)),S2323]
       
        S_3131=[ 1,(4-5*v)/(15*(1-v)),0]

        K_r=0
        
        L_r=0
        P_r=0
        M_r=0
        N_r=0
        
        k_0=E0/(3*(1-2*nu0))
        mu_0=E0/(2*(1+v))
        
        K_0=k_0+(1/3)*(mu_0)
        
        M_0=mu_0
        P_0=mu_0
        L_0=k_0-(2/3)*(mu_0)
        
        N_0=k_0+(4/3)*(mu_0)
        E_0=E0
        E0=E_0
        M0=M_0
        E_1=E1
        k_1=E1/(3*(1-2*nu1))
        mu_1=E1/(2*(1+nu1))
        
        K_1=k_1+(1/3)*(mu_1)
       
        M_1=mu_1
        P_1=mu_1
        L_1=k_1-(2/3)*(mu_1)
        
        N_1=k_1+(4/3)*(mu_1)
        E_1=E1
        E1=E_1
        M1=M_1
        
        

        K_r=[K_0, K_1,K_1]
        E_r=[E_0, E_1,E_1]
        L_r=[L_0 ,L_1,L_1]
        N_r=[N_0 ,N_1,N_1]
        P_r=[P_0 ,P_1,P_1]
        M_r=[M_0 ,M_1,M_1]

        c_r=np.empty(mm)
        d_r=np.empty(mm)
        e_r=np.empty(mm)
        f_r=np.empty(mm)
        g_r=np.empty(mm)
        h_r=np.empty(mm)
        d_CLA=np.empty(mm)
        e_CLA=np.empty(mm)
        f_CLA=np.empty(mm)
        g_CLA=np.empty(mm)
        h_CLA=np.empty(mm)
        
        c_CLA=np.empty(mm)
        term1=np.empty(mm)
        term2=np.empty(mm)
        term3=np.empty(mm)
        term4=np.empty(mm)
        
        c_CAI=np.empty(mm)
        d_CAI=np.empty(mm)
        e_CAI=np.empty(mm)
        f_CAI=np.empty(mm)
        g_CAI=np.empty(mm)
        h_CAI=np.empty(mm)
        C_r=[1-den1-den2,den1,den2]
        l_r=np.empty(mm)
       # mm=2
        
        for r in np.arange(0,mm):
            c_r[r]=1+(2*(K_r[r]-K_0)/E0)*((1-v)*(S_2222[r]+S_2233[r])-2*v*S_2211[r])+(2*(L_r[r]-L_r[0])/E0)*((1-v)*S_1122[r]-v*S_1111[r])#((2*(K_r[r]-K_r[0])/E_r[0])*(((1-v)*(S_2222[r]+S_2233[r]))-2*v*S_2211[r]))+2*((L_r[r]-L_r[0])/E_r[0])*((1-v)*(S_1122[r])-v*S_1111[r])#+c_r #Checked
            c_r[1]=1+(1/((15*(1-v))*M_0))*(2*(3-5*v)*(K_r[1]-K_0)-(L_r[1]-L_0))
            d_r[r]=(1+((N_r[r]-N_r[0])/E_r[0])*(S_1111[r]-2*v*S_1122[r])+2*((L_r[r]-L_r[0])/E_r[0])*((1-v)*(S_1122[r])-v*S_1111[r]))#+d_r #Checked
            e_r[r]=(1+2*((M_r[r]-M_r[0])/M_r[0])*S_2323[r])#+e_r
            
            f_r[r]=(1+2*S_1212[r]*((P_r[r]-P_r[0])/M_r[0]))#+f_r
           # print(f_r[r])
            g_r[r]=(2*((K_r[r]-K_r[0])/E_r[0])*((1-v)*S_1122[r]-v*S_1111[r])+1*((L_r[r]-L_r[0])/E_r[0])*((S_1111[r])-2*v*S_1122[r]))#+g_r   #checked
            h_r[r]=((N_r[r]-N_r[0])/E_r[0])*(S_2211[r]-v*(S_2222[r]+S_2233[r]))+((L_r[r]-L_r[0])/E_r[0])*(((1-v)*(S_2233[r]+S_2222[r]))-2*v*S_2211[r])#+h_r #Checked
            h_r[1]=(1/(30*(1-v)*M_0))*(-(N_r[1]-N_0)+2*(3-5*v)*(L_r[r]-L_0))
            l_r[r]=(c_r[r]*d_r[r]-2*g_r[r]*h_r[r])
            d_CLA[r]=((C_r[r]*(N_r[r]*c_r[r]-2*L_r[r]*h_r[r]))/l_r[r])#d_CLA checked
            e_CLA[r]=(2*C_r[r]*M_r[r])/e_r[r]
            f_CLA[r]=(2*C_r[r]*P_r[r]/f_r[r])
            g_CLA[r]=(C_r[r]*(L_r[r]*d_r[r]-N_r[r]*g_r[r])/l_r[r])#g_cla checked
            h_CLA[r]=(C_r[r]*(L_r[r]*c_r[r]-2*K_r[r]*h_r[r])/l_r[r]) #h_cla checked
            c_CLA[r]=((2*C_r[r]*(K_r[r]*d_r[r]-L_r[r]*g_r[r]))/l_r[r])   #c_CLA checked.
            term1[r]=(C_r[r]*d_r[r]/l_r[r])
            term2[r]=(C_r[r]*c_r[r]/l_r[r])
            term3[r]=(C_r[r]*g_r[r]/l_r[r])
            term4[r]=(C_r[r]*h_r[r]/l_r[r])
            c_CAI[r]=((C_r[r]*c_r[r]/l_r[r]))
            d_CAI[r]=((C_r[r]*d_r[r]/l_r[r]))
            e_CAI[r]=((C_r[r]/e_r[r]))
            
            f_CAI[r]=(C_r[r]/f_r[r])
            g_CAI[r]=((C_r[r]*g_r[r]/l_r[r]))
            h_CAI[r]=((C_r[r]*h_r[r]/l_r[r]))

        l_CAI=np.sum(term1)*np.sum(term2)-2*np.sum(term3)*np.sum(term4)#[r]*term2[r]-2*term3[r]*term4[r]
        d_CLA2=np.sum(d_CLA)
        e_CLA2=np.sum(e_CLA)
        f_CLA2=np.sum(f_CLA)
        g_CLA2=np.sum(g_CLA)
        h_CLA2=np.sum(h_CLA)
        c_CLA2=np.sum(c_CLA)
        
        c_CAI2=np.sum(c_CAI)/l_CAI
        d_CAI2=np.sum(d_CAI)/l_CAI
        e_CAI2=1/np.sum((e_CAI))
        
        f_CAI2=1/np.sum((f_CAI))
     #   print(f_CAI2)
        g_CAI2=np.sum(g_CAI)/l_CAI
        h_CAI2=np.sum(h_CAI)/l_CAI

        
        k_2=c_CLA2*c_CAI2+2*h_CLA2*g_CAI2 

        l=g_CLA2*c_CAI2+d_CLA2*g_CAI2


        l_=h_CLA2*d_CAI2+c_CLA2*h_CAI2

        
        
        n=d_CLA2*d_CAI2+2*g_CLA2*h_CAI2
        m_2=e_CLA2*e_CAI2
        p_2=f_CLA2*f_CAI2
        L=(k_2,l,l_,n,m_2,p_2)
        E11=n-l_**2/(k_2/2)
        mu23=(m_2/2)#/M_0
        mu12=p_2/2
        
        K23=k_2/2

        v12=l_/k_2
        E221=(1/4)*((1/(mu23)+1/(K23)))+(v12**2)/E11
        E22=(1/(E221))
        m=(m_2)

        nu12 = (E11/E22 - E11*(1./mu23 + 1./K23)/4.)**.5
        nu23 = E22/2./mu23 - 1.
        S11 = 1/E22
        S33 = 1/E11
        S12 = -nu23/E22
        S23 = -nu12/E11
        S31 = S23  
        S13 = S23
        S44 = 1/mu12
        S55 = S44
        S66 = 1/mu23
       
        C11 = S33/(2*(- 2*S13**2 + S11*S33 + S12*S33)) + 1/(2*(S11 - S12))
        C33 = -(S11 + S12)/(2*S13**2 - S33*(S11 + S12))
        C12 = S33/(2*(- 2*S13**2 + S11*S33 + S12*S33)) - 1/(2*(S11 - S12))
        C13 = S13/(2*S13**2 - S33*(S11 + S12))
        C44 = 1/S44
        C66 = 1/S66
        return [C11, C33, C12, C13, C44, C66]
       
    def spheroid(self,E0=1500., nu0=0.32, E1=73000., nu1=0.22, alpha=30., c=0.35, type='prolate'):
        
        # Convert input elastic modulus and poisson's ratio to lame' constant and shear modulus
        lambda0 = E0*nu0/((1.+nu0)*(1.-2.*nu0))
        mu0 = E0/(2.*(1.+nu0))
        lambda1 = E1*nu1/((1.+nu1)*(1.-2.*nu1))
        mu1 = E1/(2.*(1.+nu1))
        
        if type == 'prolate':
            g = alpha/(alpha**2. - 1.)**1.5*(alpha*(alpha**2. - 1.)**.5 - arccosh(alpha))       #prolate
        elif type == 'oblate':
            g = alpha/(1. - alpha**2.)**1.5*(arccos(alpha) - alpha*(1. - alpha**2.)**.5)        #oblate 
        
        D1 = 1. + 2.*(mu1 - mu0)/(lambda1 - lambda0)
        D2 = (lambda0 + 2.*mu0)/(lambda1 - lambda0)
        D3 = lambda0/(lambda1 - lambda0)
        
        S1111 = 1./(2.*(1. - nu0))*(1. - 2.*nu0 + (3.*alpha**2. - 1.)/(alpha**2. - 1.) - (1. - 2.*nu0 + 3.*alpha**2./(alpha**2. - 1.))*g)
        
        S2222 = 3./(8.*(1. - nu0))*alpha**2./(alpha**2. - 1.) + 1./(4.*(1. - nu0))*(1. - 2.*nu0 - 9./(4.*(alpha**2. - 1.)))*g
        S3333 = S2222
        
        S2233 = 1./(4.*(1. - nu0))*(alpha**2./(2.*(alpha**2. - 1.)) - (1. - 2.*nu0 + 3./(4.*(alpha**2 - 1.)))*g)
        S3322 = S2233
        
        S2211 = -1./(2.*(1. - nu0))*alpha**2./(alpha**2. - 1.) + 1./(4*(1. - nu0))*(3*alpha**2./(alpha**2. - 1.) - (1. - 2.*nu0))*g
        S3311 = S2211
        
        S1122 = -1./(2.*(1. - nu0))*(1. - 2.*nu0 + 1./(alpha**2. - 1.)) + 1./(2.*(1. - nu0))*(1. - 2.*nu0 + 3./(2.*(alpha**2. - 1.)))*g
        S1133 = S1122
        
        S2323 = 1./(4.*(1.-nu0))*(alpha**2./(2.*(alpha**2. - 1.)) + (1. - 2.*nu0 - 3./(4.*(alpha**2. - 1.)))*g)
        S3232 = S2323
        
        S1212 = 1./(4.*(1.-nu0))*(1. - 2.*nu0 - (alpha**2. + 1.)/(alpha**2. - 1.) - .5*(1. - 2.*nu0 - (3.*(alpha**2. + 1.))/(alpha**2. - 1.))*g)
        S1313 = S1212
        
        B1 = c*D1 + D2 + (1. - c)*(D1*S1111 + 2.*S2211)
        B2 = c + D3 + (1. - c)*(D1*S1122 + S2222 + S2233)
        B3 = c + D3 + (1. - c)*(S1111 + (1. + D1)*S2211)
        B4 = c*D1 + D2 + (1. - c)*(S1122 + D1*S2222 + S2233)
        B5 = c + D3 + (1. - c)*(S1122 + S2222 + D1*S2233)
        
        A1 = D1*(B4 + B5) - 2.*B2
        A2 = (1. + D1)*B2 - B4 - B5
        A3 = B1 - D1*B3
        A4 = (1. + D1)*B1 - 2.*B3
        A5 = (1. - D1)/(B4 - B5)
        A = 2.*B2*B3 - B1*(B4 + B5)
        
        # All the commented lines below are for stress calculation if a strain state were to be provided
        # e11s = (A1*e110 - A2*(e220 + e330))/A                               #eq 18 from ref
        # e22s = (2*A3*e110 + (A4 + A5*A)*e220 + (A4 - A5*A)*e330)/(2.*A)     #eq 19 from ref
        # e33s = (2*A3*e110 + (A4 - A5*A)*e220 + (A4 + A5*A)*e330)/(2.*A)     #eq 20 from ref
        # e12s = -e120/(mu0/(mu1 - mu0) + c + 2.*(1. - c)*S1212)              #eq 21 from ref
        # e23s = -e230/(mu0/(mu1 - mu0) + c + 2.*(1. - c)*S2323)              #eq 22 from ref
        # e13s = -e130/(mu0/(mu1 - mu0) + c + 2.*(1. - c)*S1313)              #no equation provided in text, but a clear description makes the equation obvious
        
        # #eq 6 from ref
        # e11pt = S1111*e11s + S1122*e22s + S1133*e33s
        # e22pt = S2211*e11s + S2222*e22s + S2233*e33s
        # e33pt = S3311*e11s + S3322*e22s + S3333*e33s
        # e12pt = S1212*e12s
        # e23pt = S2323*e23s
        # e13pt = S1313*e13s
        
        # #eq 9 from ref
        # e11t = -c*(e11pt - e11s)
        # e22t = -c*(e11pt - e11s)
        # e33t = -c*(e11pt - e11s)
        # e12t = -c*(e12pt - e12s)
        # e23t = -c*(e23pt - e23s)
        # e13t = -c*(e13pt - e13s)
        
        # #eq 4 from ref (average stress in the matrix)
        # C11110 = lambda0 + 2.*mu0
        # C22220 = lambda0 + 2.*mu0
        # C33330 = lambda0 + 2.*mu0
        # C11220 = lambda0
        # C22110 = C11220
        # C11330 = lambda0
        # C33110 = C11330
        # C22330 = lambda0
        # C33220 = C22330
        # C12120 = mu0
        # C23230 = mu0
        # C13130 = mu0
        # s11m = C11110*(e110 + e11t) + C11220*(e220 + e22t) + C11330*(e330 + e33t)
        # s22m = C22110*(e110 + e11t) + C22220*(e220 + e22t) + C22330*(e330 + e33t)
        # s33m = C33110*(e110 + e11t) + C33220*(e220 + e22t) + C33330*(e330 + e33t)
        # s23m = C23230*(e230 + e23t)
        # s12m = C12120*(e120 + e12t)
        # s13m = C13130*(e130 + e13t)
        
        # #eq 5 from ref (average stress in the inclusion)
        # C11111 = lambda1 + 2.*mu1
        # C22221 = lambda1 + 2.*mu1
        # C33331 = lambda1 + 2.*mu1
        # C11221 = lambda1
        # C22111 = C11221
        # C11331 = lambda1
        # C33111 = C11331
        # C22331 = lambda1
        # C33221 = C22331
        # C12121 = mu1
        # C23231 = mu1
        # C13131 = mu1
        # s11i = C11111*(e110 + e11t + e11pt) + C11221*(e220 + e22t + e22pt) + C11331*(e330 + e33t + e33pt)
        # s22i = C22111*(e110 + e11t + e11pt) + C22221*(e220 + e22t + e22pt) + C22331*(e330 + e33t + e33pt)
        # s33i = C33111*(e110 + e11t + e11pt) + C33221*(e220 + e22t + e22pt) + C33331*(e330 + e33t + e33pt)
        # s23i = C23231*(e230 + e23t + e23pt)
        # s12i = C12121*(e120 + e12t + e12pt)   
        # s13i = C13131*(e130 + e13t + e13pt)
        
        # Some elastic constants of the composite
        E11 = E0/(1. + c*(A1 + 2.*nu0*A2)/A)                                            # eq25 from ref
        E22 = E0/(1. + c*(-2.*nu0*A3 + (1. - nu0)*A4 + (1. + nu0)*A5*A)/(2.*A))         # eq28 from ref
        mu12 = mu0*(1. + c/(mu0/(mu1-mu0) + 2.*(1.-c)*S1212))
        mu23 = mu0*(1. + c/(mu0/(mu1-mu0) + 2.*(1.-c)*S2323))
        f = lambda K23: (lambda0+mu0)*((1.+nu0)*(1.-2.*nu0)/(1.-nu0*(1.+2.*(E11/E22 - E11*(1./mu23 + 1./K23)/4.)**.5) + c*(2.*((E11/E22 - E11*(1./mu23 + 1./K23)/4.)**.5-nu0)*A3 + (1.-nu0*(1.+2.*(E11/E22 - E11*(1./mu23 + 1./K23)/4.)**.5))*A4)/A)) - K23
        # K23 = scipy.optimize.bisect(f, E0*1.5, 2.*min((E0,E1)))
        K23 = scipy.optimize.newton(f, (1.-c)*(lambda0+mu0) + c*(lambda1+mu1))
        nu12 = (E11/E22 - E11*(1./mu23 + 1./K23)/4.)**.5
        nu23 = E22/2./mu23 - 1.
       
        m_2=mu23*2
        p_2=mu12*2
        k_2=K23*2
        E221=1/E22
        
        v12=E11/E22-(E11/4)*((1/mu23)+(1/K23))#(((1/4)*((1/(mu23)+1/(K23)))-E221)*E11)**0.5
        l_=v12*k_2
        n=E11+l_**2/(k_2/2)
        S11 = 1/E22
        S33 = 1/E11
        S12 = -nu23/E22
        S23 = -nu12/E11
        S31 = S23  
        S13 = S23
        S44 = 1/mu12
        S55 = S44
        S66 = 1/mu23
        
        C11 = S33/(2*(- 2*S13**2 + S11*S33 + S12*S33)) + 1/(2*(S11 - S12))
        C33 = -(S11 + S12)/(2*S13**2 - S33*(S11 + S12))
        C12 = S33/(2*(- 2*S13**2 + S11*S33 + S12*S33)) - 1/(2*(S11 - S12))
        C13 = S13/(2*S13**2 - S33*(S11 + S12))
        C44 = 1/S44
        C66 = 1/S66
        
        return [C11, C33, C12, C13, C44, C66,k_2,l_,n,m_2,p_2]        
       
       
    def getangle(self,Q):
        P=0.5
        #Q=11
        g_den=lambda thetaI: np.sin(thetaI/1)**(2*P-1)*np.cos(thetaI/1)**(2*Q-1)#*2
       
        gg_den=integrate.quad(g_den,0,np.pi/2)  # Limit here needs to be consistent with that in the triple integral
        g=lambda theta:((((np.sin(theta/1))**(2*P-1))*((np.cos(theta/1))**(2*Q-1)))/(gg_den[0]))
    
        N=10000
        a_n = 0
        b_n = 90
        for n in range(1,N+1):
            m_n = (a_n + b_n)/2
            f_m_n = x=(integrate.quad(g,0,m_n*3.14/180))[0]#f(m_n)
            if f_m_n>0.9:
                a_n = a_n
                b_n = m_n
            elif f_m_n< 0.9:
                a_n = m_n
                b_n = b_n
            elif f_m_n-0.9 <0.001:
               # print("Found exact solution.")
                return m_n,n
            else:
               # print("Bisection method fails.")
                return None
        return (a_n + b_n)/2, n
    #%%
    def getQ(self,theta):
        N=10000
        a_n = 0.5
        b_n = 100
        for n in range(1,N+1):
            m_n = (a_n + b_n)/2
            f_m_n = self.getangle(m_n)[0]
          #  print(f_m_n)
          #  b_n = m_n
            if abs(f_m_n-theta )<0.1:
                #print("Found exact solution.")
                return m_n
              #  break
            if f_m_n>theta:
                a_n = m_n
                b_n = b_n
            elif f_m_n< theta:
                a_n = a_n
                b_n=m_n
            else:
                return None
        
         #   return m_n,n
    
        return (a_n + b_n)/2
       
       
      