# From Tandon, G.P. and Weng, G.J., 1984. The effect of aspect ratio of inclusions on the elastic properties of unidirectionally aligned composites. Polymer composites, 5(4), pp.327-333.

from numpy import arccosh
from numpy import arccos
from numpy import linspace
from numpy import pi
from numpy import logspace
from numpy import gradient
import scipy.optimize

#############################################################################################################################          
def moduli2orthocompliance(M):
    # Order of entries in M should go as: #     E1    E2    E3    G23   G13   G12   v21   v31   v32   v12   v13   v23
                                         #       0     1     2     3     4     5     6     7     8     9     10    11
    
    S = np.empty(9)
    S[0] = 1./M[0]
    S[1] = 1./M[1]
    S[2] = 1./M[2]
    S[3] = -M[6]/M[1]
    S[4] = -M[7]/M[2]
    S[5] = -M[8]/M[2]
    S[6] = 1./M[3]
    S[7] = 1./M[4]
    S[8] = 1./M[5]
    
    
    return(S)
    # Order of entries in S go as: S11, S22, S33, S12, S13, S23, S44, S55, S66
    # indices:                      0    1    2     3    4    5    6    7   8  
#############################################################################################################################
def spheroid(E0=1500., nu0=0.32, E1=73000., nu1=0.22, alpha=30., c=0.35, type='prolate'):
    
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
    print('K0')
    print(lambda0+mu0)
    print('K1')
    print(lambda1+mu1)
    print('K23')
    print(K23)
    print('nu12')
    print(nu12)
    print('E11')
    print(E11)
    print('E22')
    print(E22)
    print('mu12')
    print(mu12)
    print('mu23')
    print(mu23)
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
     