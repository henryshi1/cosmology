import math
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
#                               Function dictionary
# ----------------------------------------------------------------------------

# Calculate E(z)
# This function is implemented badly for density_m = 0.05
def E(z_i, density_m, density_l):
    density_k = 0
    E = np.sqrt(density_m * (1.0 + z_i)**3 + density_k * (1.0 + z_i)**2 + density_l)
    return E;

# Calculate the integral D_C/D_H, store it in an array
def DCDH_int(z_limit, density_m, density_l):
    f = lambda z_i: 1.0 / E(z_i, density_m, density_l)
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# for some reason, these functions try to make DMDH a scalar, even though DMDH is an array
# proper motion distance for density_l == 0
def prop_motion_0(z, DMDH, density_m, density_l):
    # create the output array DMDH
    DMDH = np.copy(z);

    # find density_k = 1 - density_m - density_l
    density_k = 1.0 - density_m - density_l
    print(r"$\Omega_k =$ " + str(density_k))
    print(r"$\Omega_\Mu =$ " + str(density_m))
    print(r"$\Omega_\Lambda =$ " + str(density_l))

    # Calculate D_M/D_H for different universe geometries
    # open universe
    if (density_k > 0):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(density_k) * np.sinh( np.sqrt(density_k) * DCDH_int(z[i],density_m,density_l) )
    # closed universe
    elif (density_k < 0):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(abs(density_k)) * np.sin( np.sqrt(abs(density_k)) * DCDH_int(z[i],density_m,density_l) )
    # flat universe
    else:
        DMDH = 2.0*( 2.0 - density_m*(1.0-z) - (2-density_m)*np.sqrt(1.0+density_m*z) ) / ( density_m*density_m*(1.0+z) );

    return DMDH

# proper motion distance for density_l != 0
def prop_motion(z, DMDH, density_m, density_l):
    DMDH = np.copy(z);
    
    # find density_k = 1 - density_m - density_l
    density_k = 1.0 - density_m - density_l
    print(density_k)
    
    # Calculate D_M/D_H for different universe geometries
    if (density_k > 0):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(density_k) * np.sinh( np.sqrt(density_k) * DCDH_int(z[i],density_m,density_l) )
    # closed universe
    elif (density_k < 0):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(abs(density_k)) * np.sin( np.sqrt(abs(density_k)) * DCDH_int(z[i],density_m,density_l) )
    # flat universe
    else:
        for i in range(len(z)):
            DMDH[i] = DCDH_int(z[i], density_m, density_l)
    return DMDH

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
def dist_mod(z, D_H, density_m, density_l):
    # initialize the distance modulus array
    mu = np.copy(z)

    # convert D_H from Mpc to pc
    D_H = D_H * 1.0e6

    # calculate each value of mu array
    for i in range(len(z)):
        mu[i] = 5.0 * ( np.log10(1.0+z[i]) + np.log10(DCDH_int(z[i], density_m, density_l)) + np.log10(D_H/10) )
    
    return mu;

# comoving volume element
def com_vol(z, DMDH, DVC, density_m, density_l, D_H):
    DVC = np.copy(z)
    for i in range(len(z)):
        DVC[i] = DMDH[i]**2 * D_H**3 / ( E(z[i], density_m, density_l) )
    return DVC;

# ------------------------------------------------------------------------
# lookback time and age

# lookback integral
def lookback_integral(z_limit, density_m, density_l, H):
    f = lambda z: 1.0 / ((1.0+z) * E(z, density_m, density_l))
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# lookback time
def lookback_time(z, tLtH, density_m, density_l, H):
    tLtH = np.copy(z)
    for  i in range(len(z)):
        tLtH[i] = lookback_integral(z[i], density_m, density_l, H)
    return tLtH;

# age

# ------------------------------------------------------------------------

# intersection probability
def intersect_prob(z, dPdz, density_m, density_l):
    dPdz = np.copy(z)
    for i in range(len(z)):
        dPdz[i] = (1.0+z[i]) * (1.0+z[i]) /  E(z[i], density_m, density_l)
    return dPdz;

# ------------------- Custom functions that change global variables --------------------

def calculate_DMDH():
    if (density_l == 0):
        DMDH = prop_motion_0(z, DMDH, density_m, density_l)
    else:
        DMDH = prop_motion(z, DMDH, density_m, density_l)
    return

# --------------------------------------------------------------------------------------
#                                   Variable dictionary
# --------------------------------------------------------------------------------------

# Constants
density_k = 0                   # critical density
c = 2.9979e5                    # speed of light (km/s)
density_m = None                # mass density parameter
density_l = None                # energy density

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
param_text = ""                 # label for the plots

# Output arrays
z = None;                       # Array of input values of z
DMDH = np.array([]);            # Array of output values of D_M/D_H
DADH = np.array([]);            # Array of output values of D_A/D_H
DLDH = np.array([]);            # Array of output values of D_L/D_H
dmod = np.array([]);            # distance modulus array
DVC = np.array([]);             # comoving volume array
tLtH = np.array([]);            # Array of lookback times
dPdz = np.array([]);            # dimensionless intersection probability

# arrays of density parameters
omega_m = np.array([1.0, 0.05, 0.2])
omega_l = np.array([0.0, 0.0, 0.8])

# ---------------------------------------------------------------------------------------------
#                                           User input
# ---------------------------------------------------------------------------------------------

# Take in densities
print("\u03A9_k = 0.0")

# Take in Hubble constant
H = 70#float(input("Enter Hubble constant (in km/s/Mpc): "))
D_H = c / H                                         # Hubble distance

# Parameters for plot, including the range of redshift values and y-axis
zmin = 0
zmax = 5
ymin = 0
ymax = None

# Define the array of z-values (redshift)
z = np.linspace(zmin, zmax, 100)
print("z values:");
print(z);

# ----------------------------------------------------------------------------------------------
#                     Choose a option and it will calculate and plot for you                    
# ----------------------------------------------------------------------------------------------

# Ask user for the desired plot
plot_type = (input("Enter the number plot you want with respect to z, as an integer"
                    + "\n" + "1. proper motion distance D_M/D_H"
                    + "\n" + "2. angular diameter distance D_A/D_H"
                    + "\n" + "3. luminosity distance D_L/D_H"
                    + "\n" + "4. distance modulus DM + 5log(h)"
                    + "\n" + "5. comoving volume element [1/D_H]^3 dV_C/dz/d\u03A9"
                    + "\n" + "6. lookback time t_L/t_H and age t/t_H"
                    + "\n" + "7. dimensionless intersection probability dP/dz"
                    + "\n"))
option = int(plot_type)

# Ask user for maximum y-axis value
ymax = float(input("ymax: "));

#Display plot - universal aspects of every plot
plt.xlabel("redshift z");
plt.legend(loc='lower right');
plt.axis([zmin, zmax, ymin, ymax]);

# 1. proper motion distance
if (option == 1):
    print("Proper Motion Distance vs Redshift");
    plt.title("Proper Motion Distance vs Redshift");
    plt.ylabel(r"proper motion distance $D_M/D_H$");

    # Note: I want to create a new version of DMDH with each loop iteration
    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        if (density_l == 0):
            DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        else:
            DMDH = prop_motion(z, DMDH, density_m, density_l)
        print("DMDH")
        print(DMDH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, DMDH, label=param_text)

# 2. angular diameter distance
elif (option == 2):
    print("Angular Diameter Distance vs Redshift");
    plt.title("Angular Diameter Distance vs Redshift");
    plt.ylabel(r"angular diameter distance $D_A/D_H$");

    # Note: I want to create a new version of DADH with each loop iteration
    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        if (density_l == 0):
            DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        else:
            DMDH = prop_motion(z, DMDH, density_m, density_l)
        DADH = dist_angdiam(z, DMDH, DADH, density_m, density_l)
        print("DADH")
        print(DADH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, DADH, label=param_text)
    
# 3. luminosity distance
elif (option == 3):
    print("Luminosity Distance vs Redshift");
    plt.title("Luminosity Distance vs Redshift");
    plt.ylabel(r"luminosity distance $D_L/D_H$");

    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        if (density_l == 0):
            DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        else:
            DMDH = prop_motion(z, DMDH, density_m, density_l)
        DLDH = dist_lum(z, DMDH, DLDH, density_m, density_l)
        print("DLDH")
        print(DLDH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, DLDH, label=param_text)

# 4. distance modulus
elif (option == 4):
    print("Distance Modulus vs Redshift");
    plt.title("Distance Modulus vs Redshift");
    plt.ylabel("distance modulus DM + 5 log h (mag)");

    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        if (density_l == 0):
            DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        else:
            DMDH = prop_motion(z, DMDH, density_m, density_l)
        mu = dist_mod(z, D_H, density_m, density_l)
        print("i = " + str(i))
        print("distance molulus DM + 5 log h")
        print("length of mu: " + str(len(mu)))
        print(mu)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, DLDH, label=param_text)

# 5. comoving volume element
elif (option == 5):        
    print("Comoving Volume Element vs Redshift");
    plt.title("Comoving Volume Element vs Redshift");
    plt.ylabel(r"comoving volume element $[1/D_H]^3$" + r" $dV_C/dz/d\Omega$");

    # generates the same plot regardless of the parameter you give it
    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        if (density_l == 0):
            DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        else:
            DMDH = prop_motion(z, DMDH, density_m, density_l)
        DVC = com_vol(z, DMDH, DVC, density_m, density_l, D_H)
        DVC = DVC / D_H**3
        print("DVC")
        print(DVC)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, DVC, label=param_text)

# 6. Lookback time and age vs redshift
elif (option == 6):
    print("Lookback Time and Age vs Redshift");
    plt.title("Lookback Time and Age vs Redshift");
    plt.ylabel(r"lookback time $t_L/t_H$" + r" and age $t/t_H$");

    # lookback time vs redshift
    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        #if (density_l == 0):
        #    DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        #else:
        #    DMDH = prop_motion(z, DMDH, density_m, density_l)
        tLtH = lookback_time(z, tLtH, density_m, density_l, H)
        print("tLtH")
        print(tLtH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, tLtH, label=param_text)
        
    # age of universe vs redshift

# 7. intersection probability vs redshift
elif (option == 7):        
    print("Dimensionless Intersection Probability vs Redshift");
    plt.title("Proper Motion Distance vs Redshift");
    plt.ylabel(r"proper motion distance $\Omega_M/\Omega_H$");
    for i in range(len(omega_m)):
        density_m = omega_m[i]
        density_l = omega_l[i]
        print(density_m, density_l)
        #if (density_l == 0):
        #    DMDH = prop_motion_0(z, DMDH, density_m, density_l)
        #else:
        #    DMDH = prop_motion(z, DMDH, density_m, density_l)
        dPdz = intersect_prob(z, dPdz, density_m, density_l)
        print("dPdz")
        print(dPdz)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(density_m) + ", " + str(density_l) + ")$"
        plt.plot(z, dPdz, label=param_text)
    
else:
    print("Error: Please select an integer in the range 1-7.");

plt.legend()
plt.show()
