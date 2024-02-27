
import numpy as np
#C11, C22, C33, C12, C13, C23, C44, C55, C66 
def StoM(C):
    #C=[40.68,625.7,40.68,12.4,39.2,12.4,2.44,1.36,2.44]
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
    
    
    E = [1./S[i] for i in range(3)]   
    G = [1./S[i] for i in range(6,9)]
    nu21 = -E[1]*S[3]
    nu31 = -E[2]*S[4]
    nu32 = -E[2]*S[5]
    nu12 = -E[0]*S[3]
    nu13 = -E[0]*S[4]
    nu23 = -E[1]*S[5]
    print([E[0], E[1], E[2], G[0], G[1], G[2], nu21, nu31, nu32, nu12, nu13, nu23])