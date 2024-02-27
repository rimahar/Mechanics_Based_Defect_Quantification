# Import packages
import math as math
import numpy as np
import Two_Inclusions as md
from numpy import arccosh
from numpy import arccos
from numpy import linspace
from numpy import pi
from numpy import logspace
from numpy import gradient
import scipy.optimize
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.integrate import tplquad
import matplotlib.pyplot as plt
import eshelby_Spheroid
def ODF(E0,nu0,E1,nu1,c_volume,Q,numberOfinclusions,alpha1):
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
    L_ijkl4=np.zeros(81)
    L_ijkl4=L_ijkl3.reshape(3,3,3,3)
    L_constant1=np.zeros(81)
    L_constant1=L_constant1.reshape(3,3,3,3)
    L_constant2=np.zeros(81)
    L_constant2=L_constant2.reshape(3,3,3,3)
    
    #Defining L_pqrs
    L_pqrs=np.zeros(81)
    L_pqrs=L_pqrs.reshape(3,3,3,3)
    L_matrix = np.zeros((3,3,3,3))
    L_pqrs_phase2=np.zeros(81)
    L_pqrs_phase2=L_pqrs_phase2.reshape(3,3,3,3)
   
    
    Cmatrix = eshelby_Spheroid.spheroid(E0, nu0, E1, nu1, alpha=100., c=0.0, type='prolate')
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
    
    
    Cvoids=0
    Ccracks=1-Cvoids
    C=md.model(alpha1,Cvoids*c_volume,Ccracks*c_volume,numberOfinclusions,E0,nu0,E1,nu1)#C[6],C[7],C[8],C[9],C[10])

    ## Enter L
    # Unique entries in L
    L_pqrs[2][2][2][2] = C[1]
    L_pqrs[2][2][0][0] = C[3]
    L_pqrs[0][0][0][0] = C[0]
    L_pqrs[0][0][1][1] = C[2]
    L_pqrs[1][1][2][2] = C[3]
    L_pqrs[1][1][1][1] = C[0]
    L_pqrs[0][1][0][1] = C[5]#2
    L_pqrs[2][1][2][1] = C[4]#2
    L_pqrs[2][0][2][0] = C[4]#/2
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
    
    
    
    
    
    L_constant1[2][2][2][2] = 21853.93
    L_constant1[2][2][0][0] = 1681.11
    L_constant1[0][0][0][0] = 3701.88
    L_constant1[0][0][1][1] = 1633.87
    L_constant1[1][1][2][2] = 1681.11
    L_constant1[1][1][1][1] = 3701.88
    L_constant1[0][1][0][1] = 1034#2
    L_constant1[2][1][2][1] = 1288.58#2
    L_constant1[2][0][2][0] = 1288.58#/2
    # Symmetry of off diagonal normal terms
    L_constant1[0][0][2][2]=L_constant1[2][2][0][0]
    L_constant1[2][2][1][1]=L_constant1[1][1][2][2]
    L_constant1[1][1][0][0]=L_constant1[0][0][1][1]
    # Symmetry of shear terms
    L_constant1[0][2][0][2]=L_constant1[2][0][2][0]
    L_constant1[0][2][2][0]=L_constant1[2][0][2][0]
    L_constant1[2][0][0][2]=L_constant1[2][0][2][0]
    L_constant1[1][0][1][0]=L_constant1[0][1][0][1]
    L_constant1[1][0][0][1]=L_constant1[0][1][0][1]
    L_constant1[0][1][1][0]=L_constant1[0][1][0][1]
    L_constant1[1][2][1][2]=L_constant1[2][1][2][1]
    L_constant1[1][2][2][1]=L_constant1[2][1][2][1]
    L_constant1[2][1][1][2]=L_constant1[2][1][2][1]
    
    
    L_constant2[2][2][2][2] = 14054.498789826492
    L_constant2[2][2][0][0] = 1129.5204583682303
    L_constant2[0][0][0][0] = 2407.91338199813
    L_constant2[0][0][1][1] = 1061.99371575535
    L_constant2[1][1][2][2] = 1129.5204583682303
    L_constant2[1][1][1][1] = 2407.91338199813
    L_constant2[0][1][0][1] = 672.9548331213903#2
    L_constant2[2][1][2][1] = 873.9404524444342#2
    L_constant2[2][0][2][0] = 873.9404524444342#/2
    # Symmetry of off diagonal normal terms
    L_constant2[0][0][2][2]=L_constant2[2][2][0][0]
    L_constant2[2][2][1][1]=L_constant2[1][1][2][2]
    L_constant2[1][1][0][0]=L_constant2[0][0][1][1]
    # Symmetry of shear terms
    L_constant2[0][2][0][2]=L_constant2[2][0][2][0]
    L_constant2[0][2][2][0]=L_constant2[2][0][2][0]
    L_constant2[2][0][0][2]=L_constant2[2][0][2][0]
    L_constant2[1][0][1][0]=L_constant2[0][1][0][1]
    L_constant2[1][0][0][1]=L_constant2[0][1][0][1]
    L_constant2[0][1][1][0]=L_constant2[0][1][0][1]
    L_constant2[1][2][1][2]=L_constant2[2][1][2][1]
    L_constant2[1][2][2][1]=L_constant2[2][1][2][1]
    L_constant2[2][1][1][2]=L_constant2[2][1][2][1]
    
    LPQRS=[[L_pqrs[0][0][0][0],L_pqrs[0][0][1][1],L_pqrs[0][0][2][2],0,0,0],[L_pqrs[1][1][0][0],L_pqrs[1][1][1][1],L_pqrs[1][1][2][2],0,0,0],[L_pqrs[2][2][0][0],L_pqrs[2][2][1][1],L_pqrs[2][2][2][2],0,0,0],[0,0,0,L_pqrs[1][2][1][2],0,0],[0,0,0,0,L_pqrs[0][2][0][2],0],[0,0,0,0,0,L_pqrs[0][1][0][1]]]

    L_pqrs=L_pqrs-L_matrix



    counter=0
   
    g2=lambda theta:((((np.sin(theta/2))**(2*P-1))*((np.cos(theta/2))**(2*Q-1)))/(gg_den[0]))#*2#/(3.14/90)
    C=lambda beta,phi, theta: (g2(theta)*np.sin(theta))
    factor=integrate.tplquad(C, 0, np.pi/2, lambda theta: 0, lambda theta: 2*np.pi, lambda theta, phi: 0, lambda theta, phi: 2*np.pi)[0]#*4
    #factor=integrate.dblquad(C, 0, np.pi/2, lambda phi: 0, lambda phi: 2*np.pi)#, lambda theta, phi: 0, lambda theta, phi: 2*np.pi)[0]
    print('factor')
    print(factor)
    
    print(gg_den[0])
    c1=0.7
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
                                      
                                        L_ijkl3[i][j][k][l]=1*L_matrix[i][j][k][l]+1*(L_ijkl2[i][j][k][l])#+L_ijkl2_phase2[i][j][k][l]
                                     
                                        counter += 1
                                        
                        print(str(i)+str(j)+str(k)+str(l))

                        print(L_ijkl3[i][j][k][l])
                        
                                        

    
    return([L_ijkl3[0][0][0][0],L_ijkl3[0][0][0][0],L_ijkl3[2][2][2][2],L_ijkl3[0][0][1][1],L_ijkl3[0][0][2][2],L_ijkl3[0][0][2][2],L_ijkl3[1][2][1][2],L_ijkl3[1][2][1][2],L_ijkl3[0][1][0][1]])

