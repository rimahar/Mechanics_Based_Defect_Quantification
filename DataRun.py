import Orientation_Averaging as mtodf
import StifftoMod as SM
import numpy as np
import csv
def odf(E0,nu0,E1,nu1,c1,r1,r2,numberOfinclusions):
    jj=0
    #Set Values of Q
    Qall=[0.15,0.3,0.36,0.43,0.48,0.49,0.53,0.58,0.65,0.7,0.75,0.8,0.88,0.95,1,1.1,1.2,1.3,1.4,1.5,1.65,1.8,1.9,2.1,2.3 ,2.5  ,2.8 ,3.1,3.5 ,4,4.5 ,5 ,5.7,6.7 ,8,9,11,14 ,18,25,30 ,45 ,75,125,275,1100]#[0.5,2,5,8,10,12,15,18,22,25]
    Qall=np.flip(Qall)
    Cvoids=0 
    for Q in Qall:
            for alpha1 in np.logspace(3,1):
           
                C=mtodf.ODF(E0,nu0,E1,nu1,c1,Q,numberOfinclusions,alpha1)    
                    
                f1 = -C[1]*C[4]**2. + 2.*C[3]*C[4]*C[5] - C[0]*C[5]**2. + C[2]*(C[0]*C[1]-C[3]**2.)    
                f2 = (C[4]*C[5]-C[3]*C[2])/f1    
                f3 = (C[3]*C[5]-C[4]*C[1])/f1        
                f4 = (C[3]*C[4]-C[0]*C[5])/f1        
                        
                # Set values        
                # Order of entries in S go as: S11, S22, S33, S12, S13, S23, S44, S55, S66        
                S = np.empty(9)        
                S[0] = (C[1]*C[2]-C[5]**2.)/f1        
                S[1] = (C[0]*C[2]-C[4]**2.)/f1        
                S[2] = (C[0]*C[1]-C[3]**2.)/f1        
                S[3] = f2        
                S[4] = f3        
                S[5] = f4        
                S[6] = 1./C[6]        
                S[7] = 1./C[7]        
                S[8] = 1./C[8]        
                    
                E11=1/S[2]    
                E22=1/S[0]   #44,55,66 
                print(E11)    
                print(E22)    
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
                C00.append(C[0])    
                C11.append(C[1])    
                C22.append(C[2])    
                C33.append(C[3])    
                C44.append(C[4])    
                C55.append(C[5])    
                C66.append(C[6])    
                C77.append(C[7])    
                C88.append(C[8])    
                volume=[]
                QQ=[]
                QQ.append(Q) 
                alphas=[]
                alphas.append(alpha1)
                CVOIDS=[]            
                CVOIDS.append(Cvoids)
                E111.append(float(E11))    
                  
                E222.append(float(E22))    
                  
                volume.append(float(c1))    
                    
                SM.StoM(C)    
                    
                outputData=np.transpose([QQ,alphas,volume,E111,E222,C00,C11,C22,C33,C44,C55,C66,C77,C88])    
                if jj==0:    
                    ww='w'    
                else:    
                    ww='a'    
                jj+=1    
                outputFile = 'Data_'+str(c1)+'_'+str(Q)+'_'+str(round(alpha1,1))+' alpha.csv'    
                with open(outputFile, mode=ww) as csvfile:    
                    output = csv.writer(csvfile, delimiter=',', quotechar='|')    
                    for row in outputData:    
                        output.writerow(row)    
     
                outputFile = 'Data_'+str(c1)+'_all.csv'
                with open(outputFile, mode=ww) as csvfile:
                    output = csv.writer(csvfile, delimiter=',', quotechar='|')
                    for row in outputData:
                        output.writerow(row)

                  
    return QQ,CVOIDS,volume,E111,E222,C00,C11,C22,C33,C44,C55,C66,C77,C88 
        
        
