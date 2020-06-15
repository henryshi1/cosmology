import math
import numpy
import scipy.integrate as integrate

# Constant dictionary
density_m = 0.3                 # mass density
density_l = 0.7                 # energy density
density_k = 0                   # critical density
c = 2.9979e5                    # speed of light (km/s)

# Variable dictionary
H = None                        # Hubble constant (km/s/Mpc)
z = None                        # redshift
D_C = None                      # comoving distance, line-of-sight (m)
D_M = None                      # comoving distance, transverse (m)
D_A = None                      # angular diameter distance (m)
D_L = None                      # luminosity distance (m)
V_C = None                      # comoving volume (m^3)
D_H = None                      # Hubble Distance (Mpc)

# Take in Hubble constant and redshift
H = float(input("Enter Hubble constant (in km/s/Mpc): "))
z = float(input("Enter redshift: "))
# Calculate Hubble Distance (Mpc) from Hubble Constant
D_H = c / H


# comoving distance, line-of-sight
# I am having trouble with this integral -H.S.
f = lambda z: 1.0 / math.sqrt(density_m*(1.0+z)**3 + density_k*(1.0+z)**2 + density_l)
i = integrate.quad(f, 0, z)
D_C = D_H * i[0]
print("Comoving distance, line-of-sight (Mpc): " + str(D_C))


# comoving distance, transverse
if (density_l == 0):
    D_M = D_H * 2.0 * (2.0 - density_m*(1.0-z) - (2.0-density_m)*math.sqrt(1.0+density_m*z)) / (density_m**2*(1.0+z))
else:
    D_M = D_C
print("Comoving distance, transverse (Mpc): " + str(D_C))
    

# angular diameter distance
D_A = D_M / (1.0+z)
print("Angular diameter distance (Mpc): " + str(D_A))


# luminosity distance
D_L = (1.0+z)*D_M
print("Luminosity distance (Mpc): " + str(D_L))


# comoving volume
V_C = (4.0 * math.pi / 3.0) * (D_M)**3
print("Comoving volume (Mpc^3): " + str(V_C))
