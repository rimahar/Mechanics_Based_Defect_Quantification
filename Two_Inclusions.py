import numpy as np
from numpy import arccosh
from numpy import arccos
from numpy import linspace
from numpy import pi
from numpy import logspace
from numpy import gradient
import scipy.optimize
import math
def model(alpha,den1,den2,mm,E0,nu0,E1,nu1):
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
    print(f_CAI2)
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
   
   
   
   
   
   
   
  