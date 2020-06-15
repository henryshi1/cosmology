from math import *
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------

# Constant dictionary
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
array_size = None               # size of each array, which is the number of z values

# Output arrays
z = None;                       # Array of input values of z
DMDH = np.array([]);         # Array of output values of D_M/D_H
DADH = np.array([]);         # Array of output values of D_A/D_H
DLDH = np.array([]);         # Array of output values of D_L/D_H
dmod = np.array([]);         # distance modulus array
DVC = np.array([]);          # comoving volume array
tLtH = np.array([]);         # Array of lookback times
dPdz = np.array([]);         # dimensionless intersection probability

# ----------------------------------------------------------------------------

# Custom functions

# Calculate E(z)
def E(z_i, density_m, density_l):
    density_k = 0
    if (z_i == 0):
        return sqrt(density_m + density_l);
    E = sqrt(density_m * (1.0 + z_i)**3 + density_k * (1.0 + z_i)**2 + density_l)
    return E;

# Calculate the integral, store it in an array
def integral(z_limit, density_m, density_l):
    f = lambda z_i: 1.0 / E(z_i, density_m, density_l)
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# Calculate lookback integral
def lookback_integral(z_limit, density_m, density_l, H):
    f = lambda z: 1.0 / ((1.0+z) * E(z, density_m, density_l))
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# proper motion distance for density_l == 0
def prop_motion_0(z, DMDH, density_m, density_l):
    DMDH = np.copy(z);
    DMDH = 2.0*( 2.0 - density_m*(1.0-z) - (2-density_m)*sqrt(1.0+density_m*z) ) / ( density_m*density_m*(1.0+z) );
    return None;
# proper motion distance for density_l != 0
def prop_motion(z, DMDH, density_m, density_l):
    DMDH = np.copy(z);
    for i in range(len(z)):
        DMDH[i] = integral(z[i], density_m, density_l)
    return DMDH;

# angular diameter distance
def dist_angdiam(z, DMDH, DADH, density_m, density_l):
    DADH = np.copy(z);
    DADH = DMDH / (1.0+z)
    return DADH;

# luminosity distance
def dist_lum(z, DMDH, DLDH, density_m, density_l):
    DLDH = np.copy(z)
    DLDH = DMDH * (1.0+z)
    return DLDH;

# distance modulus

# comoving volume element
def com_vol(z, DMDH, DVC, density_m, density_l, D_H):
    DVC = np.copy(z)
    for i in range(len(z)):
        DVC[i] = DMDH[i]**2 * D_H**3 / ( E(z[i], density_m, density_l)*(1.0+z[i]) )
    return DVC;

# lookback time and age
def lookback_time(z, tLtH, density_m, density_l, H):
    tLtH = np.copy(z)
    for  i in range(len(z)):
        tLtH[i] = lookback_integral(z[i], density_m, density_l, H)
    return tLtH;

# intersection probability
def intersect_prob(z, dPdz, density_m, density_l):
    dPdz = np.copy(z)
    for i in range(len(z)):
        dPdz[i] = (1.0+z[i]) * (1.0+z[i]) /  E(z[i], density_m, density_l)
    return dPdz;

# ------------------------- User input -------------------------------

# Take in densities
print("\u03A9_k = 0.0")
density_m = 0.3#float(input("Enter \u03A9_\u039C: "))   # mass density
density_l = 0.7#float(input("Enter \u03A9_\u039B: "))   # energy density

# Take in Hubble constant
H = 70#float(input("Enter Hubble constant (in km/s/Mpc): "))
D_H = c / H                                         # Hubble constant

# Set parameters for plot, including the range of redshift values and y-axis
print("set parameters for plot")
zmin = 0#float(input("Enter minimum z-value: "))
zmax = 5#float(input("Enter maximum z-value: "))
ymin = 0#float(input("ymin: "));
ymax = 10#float(input("ymax: "));

# ---------------------- Initialize variables --------------------------

# arrays of density parameters
omega_m = np.array([1.0, 0.05, 0.2])
omega_l = np.array([0.0, 0.0, 0.8])

# Define the array of z-values (redshift)
z = np.linspace(zmin, zmax, 100)
print("z values");
print(z);

# -------------------------- Calculations ------------------------------

# Calculate all the arrays of proper motion distance, angular diameter distance, etc
# 1. Calculate proper motion distance and fill the DMDH array
if (density_l == 0.0):          
    DMDH = prop_motion_0(z, DMDH, density_m, density_l);
else:
    DMDH = prop_motion(z, DMDH, density_m, density_l);

print("D_M/D_H values");
print(DMDH);
    
# 2. Calculate angular diameter distance using values of DMDH array
DADH = dist_angdiam(z, DMDH, DADH, density_m, density_l);

print("D_A/D_H values")
print(DADH);

# 3. Calculate luminosity distance
DLDH = dist_lum(z, DMDH, DLDH, density_m, density_l);

print("D_L/D_H values")
print(DLDH)

# 4. Calculate distance modulus

# 5. Calculate comoving volume element
DVC = com_vol(z, DMDH, DVC, density_m, density_l, D_H);
DVC = DVC / D_H**3

print("comoving volume elements values")
print(DVC)

# 6. Calculate lookback time and age
tLtH = lookback_time(z, tLtH, density_m, density_l, H)

print("lookback time values")
print(tLtH)

# 7. Calculate dimensionless differential intersection probability
dPdz = intersect_prob(z, dPdz, density_m, density_l)

print("intersection probability values")
print(dPdz)

# -------------------------- Plotting the data ----------------------------

# Ask user for the desired plot
plot_type = (input("Enter the number plot you want with respect to z, as an integer"
                    + "\n" + "1. proper motion distance D_M/D_H"
                    + "\n" + "2. angular diameter distance D_A/D_H"
                    + "\n" + "3. luminosity distance DLA/D_H"
                    + "\n" + "4. distance modulus DM + 5log(h)"
                    + "\n" + "5. comoving volume element [1/D_H]^3 dV_C/dz/d\u03A9"
                    + "\n" + "6. lookback time t_L/t_H and age t/t_H"
                    + "\n" + "7. dimensionless intersection probability dP/dz"
                    + "\n"))
option = int(plot_type)

#Display plot
param_text = "(\u03A9_\u039C,\u03A9_\u039B)" + "=(" + str(density_m) + "," + str(density_l) + ")"
plt.xlabel("redshift z");
plt.legend(loc='bottom right');


if (option == 1):               # proper motion distance
    print("Proper Motion Distance vs Redshift");
    plt.title("Proper Motion Distance vs Redshift");
    plt.ylabel("proper motion distance D_M/D_H");
    plt.axis([zmin, zmax, ymin, ymax]);
    plt.plot(z, DMDH, '-r', label=param_text);
    plt.show()
elif (option == 2):             # angular diameter distance
    print("Angular Diameter Distance vs Redshift");
    plt.title("Angular Diameter Distance vs Redshift");
    plt.ylabel("angular diameter distance D_A/D_H");
    plt.axis([zmin, zmax, ymin, ymax]);
    plt.plot(z, DADH, '-r', label=param_text);
    plt.show()
elif (option == 3):        
    print("Luminosity Distance vs Redshift");
    plt.title("Luminosity Distance vs Redshift");
    plt.ylabel("luminosity distance D_L/D_H");
    plt.axis([zmin, zmax, ymin, ymax]);
    plt.plot(z, DLDH, '-r', label=param_text);
    plt.show()
elif (option == 4):
    print("Distance Modulus vs Redshift");
    plt.title("Distance Modulus vs Redshift");
    plt.ylabel("distance modulus DM + 5 log h (mag)");
    plt.plot(z, DMDH, '-r', label=param_text);
elif (option == 5):        
    print("Comoving Volume Element vs Redshift");
    plt.title("Comoving Volume Element vs Redshift");
    plt.ylabel("comoving volume element [1/D_H]^3 dV_C/dz/d\u03A9");
    plt.plot(z, DVC, '-r', label=param_text);
elif (option == 6):
    print("Lookback Time and Age vs Redshift");
    plt.title("Lookback Time and Age vs Redshift");
    plt.ylabel("lookback time t_L/t_H and age t/t_H");
    plt.plot(z, tLtH, '-r', label=param_text);
elif (option == 7):        
    print("Proper Motion Distance vs Redshift");
    plt.title("Proper Motion Distance vs Redshift");
    plt.ylabel("proper motion distance D_M/D_H");
    plt.plot(z, DMDH, '-r', label=param_text);
else:
    print("Error: Please select an integer in the range 1-7.");

# Display the plot
plt.show();
