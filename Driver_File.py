import DataRun as Dr
import numpy as np
import os as os

porosity=0.0
porosity=porosity/100
porosity=round(porosity,4)
jobName='Data_'+str(porosity)

try:
    os.mkdir(jobName)
except:
    os.chdir(jobName)
#Input Modulus
E0=70000
#Input Poisson
nu0=0.33
#Voids Data
E1=0  
nu1=0
#Volume fraction Range
r1=0
r2=1.01 
#Number of Inlusions       
numberOfinclusions=2 
KK=Dr.odf(E0,nu0,E1,nu1,porosity,r1,r2,numberOfinclusions+1)
outputData=np.transpose([KK[0],KK[1],KK[2],KK[3],KK[4],KK[5],KK[6],KK[7],KK[8],KK[9],KK[10],KK[11],KK[12],KK[13]])    
if jj==0:    
	ww='w'    
else:    
	ww='a'    
jj+=1    
outputFile = 'Data_'+str(porosity)+'.csv'    
with open(outputFile, mode=ww) as csvfile:    
	output = csv.writer(csvfile, delimiter=',', quotechar='|')   
	for row in outputData:    
		output.writerow(row)    
