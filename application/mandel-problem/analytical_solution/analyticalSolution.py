
G = E / (2.0 * (1+mu)) # Shear modulus [Pa]
lame = 2.0 * G * mu / (1 - 2.0 * mu) # Lame constant [Pa]

k_f = 100 * 9.869233e-13 # 100D [m2]
c_x = k_f * M * (lame + 2.0 * G)/(lame + 2.0 * G + alpha*alpha*M)

A1 = 2.0 * alpha + 2.0 * ((lame + G)/ M / alpha)
A2 = 2.0 * alpha * G / (lame + 2.0 * G)
A1_2 = A1/A2

beta_values = getBetaValue(A1_2)

def ux(F,a,lame,G,M,alpha,beta_values_,c_x,x,t):
    Mxx = lame + 2.0 * G
    Mzz = lame + 2.0 * G
    Mxz = lame
    term1 = F/a*Mxz/(Mxx*Mzz-Mxz*Mxz)
    term2 = 2F/a*()