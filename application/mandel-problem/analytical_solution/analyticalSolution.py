from functions import *

#################################################
#############Input###############################
#################################################
alpha = 1   # Biot's coefficient
M = 1.65e10 # Biot modulus [Pa]
mu = 0.4    # Poisson's ratio
E = 5.94e9  # Young's modulus [Pa]

G = E / (2.0 * (1+mu)) # Shear modulus [Pa]
lame = 2.0 * G * mu / (1 - 2.0 * mu) # Lame constant [Pa]

k_f = 100 * 9.869233e-13 # 100D [m2]
c_x = k_f * M * (lame + 2.0 * G)/(lame + 2.0 * G + alpha*alpha*M)

A1 = 2.0 * alpha + 2.0 * ((lame + G)/ M / alpha)
A2 = 2.0 * alpha * G / (lame + 2.0 * G)
A1_2 = A1/A2

def sum_beta(beta_values_list):
    sum_ = 0
    for val in beta_values_list:
        sum_ = sum_ + \
               math.sin(val) * math.cos(val)/(val - math.sin(val) * math.cos(val))
    return sum_


beta_values = list()
for i in range(50):

    beta = rn(Rn,A1_2)
    #value = math.tan(beta)/beta - A1_2
    beta_values.append(beta)

print(sum_beta(beta_values))
