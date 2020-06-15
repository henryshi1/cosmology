import math
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
#                               Function dictionary
# ----------------------------------------------------------------------------

# Calculate E(z)
def E(z_i, omega_m, omega_l):
    omega_k = 1.0 - omega_m - omega_l
    E = np.sqrt(omega_m * (1.0 + z_i)**3 + omega_k * (1.0 + z_i)**2 + omega_l)
    return E;

# Calculate the integral D_C/D_H, store it in an array
def DCDH_int(z_limit, omega_m, omega_l):
    f = lambda z_i: 1.0 / E(z_i, omega_m, omega_l)
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# proper motion distance for omega_l == 0
def prop_motion_0(z, DMDH, omega_m, omega_l):
    # create the output array DMDH
    DMDH = np.copy(z);

    # find omega_k = 1 - omega_m - omega_l
    omega_k = 1.0 - omega_m - omega_l
    print(r"$\Omega_k =$ " + str(omega_k))
    print(r"$\Omega_\Mu =$ " + str(omega_m))
    print(r"$\Omega_\Lambda =$ " + str(omega_l))

    # Calculate D_M/D_H for different universe geometries
    # open universe
    if (omega_k > 1.0e-6):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(omega_k) * np.sinh( np.sqrt(omega_k) * DCDH_int(z[i],omega_m,omega_l) )
    # closed universe
    elif (omega_k < -1.0e-6):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(abs(omega_k)) * np.sin( np.sqrt(abs(omega_k)) * DCDH_int(z[i],omega_m,omega_l) )
    # flat universe
    else:
        DMDH = 2.0*( 2.0 - omega_m*(1.0-z) - (2-omega_m)*np.sqrt(1.0+omega_m*z) ) / ( omega_m*omega_m*(1.0+z) );

    return DMDH

# proper motion distance for omega_l != 0
def prop_motion(z, DMDH, omega_m, omega_l):
    DMDH = np.copy(z);
    
    # find omega_k = 1 - omega_m - omega_l
    omega_k = 1.0 - omega_m - omega_l
    print(omega_k)
    
    # Calculate D_M/D_H for different universe geometries
    if (omega_k > 1.0e-6):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(omega_k) * np.sinh( np.sqrt(omega_k) * DCDH_int(z[i],omega_m,omega_l) )
    # closed universe
    elif (omega_k < -1.0e-6):
        for i in range(len(z)):
            DMDH[i] = 1/np.sqrt(abs(omega_k)) * np.sin( np.sqrt(abs(omega_k)) * DCDH_int(z[i],omega_m,omega_l) )
    # flat universe
    else:
        for i in range(len(z)):
            DMDH[i] = DCDH_int(z[i], omega_m, omega_l)
    return DMDH

# angular diameter distance
def dist_angdiam(z, DMDH, DADH, omega_m, omega_l):
    DADH = np.copy(z);
    DADH = DMDH / (1.0+z)
    return DADH;

# luminosity distance
def dist_lum(z, DMDH, DLDH, omega_m, omega_l):
    DLDH = np.copy(z)
    DLDH = DMDH * (1.0+z)
    return DLDH;

# distance modulus
def dist_mod(z, D_H, omega_m, omega_l, DMDH):
    # initialize the distance modulus array
    mu = np.copy(z)

    # convert D_H from Mpc to pc
    D_H = D_H * 1.0e6

    # calculate each value of mu array
    #for i in range(len(z)):
        #mu[i] = 5.0 * ( np.log10(1.0+z[i]) + np.log10(DCDH_int(z[i], omega_m, omega_l)) + np.log10(D_H/10) )
    mu = 5.0 * ( np.log10(1.0+z) + np.log10(DMDH) + np.log10(D_H/10) )
    return mu;

# comoving volume element
def com_vol(z, DMDH, DVC, omega_m, omega_l, D_H):
    DVC = np.copy(z)
    for i in range(len(z)):
        DVC[i] = DMDH[i]**2 * D_H**3 / ( E(z[i], omega_m, omega_l) )
    return DVC;

# ------------------------------------------------------------------------
# lookback time and age

# lookback integral
def lookback_integral(z_limit, omega_m, omega_l, H):
    f = lambda z: 1.0 / ((1.0+z) * E(z, omega_m, omega_l))
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# lookback time
def lookback_time(z, tLtH, omega_m, omega_l, H):
    tLtH = np.copy(z)
    for  i in range(len(z)):
        tLtH[i] = lookback_integral(z[i], omega_m, omega_l, H)
    return tLtH

# age = t(z) - t(inf)
def get_age(z, tLtH, age, omega_m, omega_l, H):
    age = np.copy(z)

    # Calculate (lookback time) / (Hubble time) at infinite redshift
    u = z / (1.0+z)
    f = lambda u: (1.0-u)**(-1) * 1.0/(E(u/(1.0-u), omega_m, omega_l))
    tInf_tH = integrate.quad(f, 0, 1)
    print(tInf_tH[0])

    # age = t(current z value) - t(infinity)
    for i in range(len(z)):
        age[i] = tInf_tH[0] - tLtH[i]

    return age

# ------------------------------------------------------------------------

# intersection probability
def intersect_prob(z, dPdz, omega_m, omega_l):
    dPdz = np.copy(z)
    for i in range(len(z)):
        dPdz[i] = (1.0+z[i]) * (1.0+z[i]) /  E(z[i], omega_m, omega_l)
    return dPdz;

# ------------------- Custom functions that change global variables --------------------

def calculate_DMDH():
    if (omega_l == 0):
        DMDH = prop_motion_0(z, DMDH, omega_m, omega_l)
    else:
        DMDH = prop_motion(z, DMDH, omega_m, omega_l)
    return

# --------------------------------------------------------------------------------------
#                                   Variable dictionary
# --------------------------------------------------------------------------------------

# Constants
omega_k = 0                     # critical density
c = 2.9979e5                    # speed of light (km/s)
omega_m = None                  # mass density parameter
omega_l = None                  # energy density

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
DMDH = np.array([])             # Array of output values of D_M/D_H
DADH = np.array([])             # Array of output values of D_A/D_H
DLDH = np.array([])             # Array of output values of D_L/D_H
dmod = np.array([])             # distance modulus array
DVC = np.array([])              # comoving volume array
tLtH = np.array([])             # Array of lookback times
age = np.array([])              # Array of ages
dPdz = np.array([])             # dimensionless intersection probability

# arrays of density parameters
omegam_array = np.array([1.0, 0.05, 0.2])
omegal_array = np.array([0.0, 0.0, 0.8])

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

#Display plot - universal aspects of every plot
plt.xlabel("redshift z");
plt.legend(loc='lower right');

# 1. proper motion distance
if (option == 1):
    print("Proper Motion Distance vs Redshift");
    plt.title("Proper Motion Distance vs Redshift");
    plt.ylabel(r"proper motion distance $D_M/D_H$");

    # Note: I want to create a new version of DMDH with each loop iteration
    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]
        if (omega_l == 0):
            DMDH = prop_motion_0(z, DMDH, omega_m, omega_l)
        else:
            DMDH = prop_motion(z, DMDH, omega_m, omega_l)
        print("DMDH")
        print(DMDH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, DMDH, label=param_text)

# 2. angular diameter distance
elif (option == 2):
    print("Angular Diameter Distance vs Redshift");
    plt.title("Angular Diameter Distance vs Redshift");
    plt.ylabel(r"angular diameter distance $D_A/D_H$");

    # Note: I want to create a new version of DADH with each loop iteration
    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]
        if (omega_l == 0):
            DMDH = prop_motion_0(z, DMDH, omega_m, omega_l)
        else:
            DMDH = prop_motion(z, DMDH, omega_m, omega_l)
        DADH = dist_angdiam(z, DMDH, DADH, omega_m, omega_l)
        print("DADH")
        print(DADH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, DADH, label=param_text)
    
# 3. luminosity distance
elif (option == 3):
    print("Luminosity Distance vs Redshift");
    plt.title("Luminosity Distance vs Redshift");
    plt.ylabel(r"luminosity distance $D_L/D_H$");

    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]
        if (omega_l == 0):
            DMDH = prop_motion_0(z, DMDH, omega_m, omega_l)
        else:
            DMDH = prop_motion(z, DMDH, omega_m, omega_l)
        DLDH = dist_lum(z, DMDH, DLDH, omega_m, omega_l)
        print("DLDH")
        print(DLDH)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, DLDH, label=param_text)

# 4. distance modulus
elif (option == 4):
    print("Distance Modulus vs Redshift");
    plt.title("Distance Modulus vs Redshift");
    plt.ylabel("distance modulus DM + 5 log h (mag)");
    ymin = 40
    ymax = 50

    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]
        if (omega_l == 0):
            DMDH = prop_motion_0(z, DMDH, omega_m, omega_l)
        else:
            DMDH = prop_motion(z, DMDH, omega_m, omega_l)
        mu = dist_mod(z, D_H, omega_m, omega_l, DMDH)
        print("i = " + str(i))
        print("distance molulus DM + 5 log h")
        print("length of mu: " + str(len(mu)))
        print(mu)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, mu, label=param_text)

# 5. comoving volume element
elif (option == 5):        
    print("Comoving Volume Element vs Redshift");
    plt.title("Comoving Volume Element vs Redshift");
    plt.ylabel(r"comoving volume element $[1/D_H]^3$" + r" $dV_C/dz/d\Omega$");

    # generates the same plot regardless of the parameter you give it
    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]
        if (omega_l == 0):
            DMDH = prop_motion_0(z, DMDH, omega_m, omega_l)
        else:
            DMDH = prop_motion(z, DMDH, omega_m, omega_l)
        DVC = com_vol(z, DMDH, DVC, omega_m, omega_l, D_H)
        DVC = DVC / D_H**3
        print("DVC")
        print(DVC)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, DVC, label=param_text)

# 6. Lookback time and age vs redshift
elif (option == 6):
    print("Lookback Time and Age vs Redshift");
    plt.title("Lookback Time and Age vs Redshift");
    plt.ylabel(r"lookback time $t_L/t_H$" + r" and age $t/t_H$");

    # lookback time and age vs redshift
    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]


        # lookback time
        tLtH = lookback_time(z, tLtH, omega_m, omega_l, H)
        print("tLtH")
        print(tLtH)
        param_text = r"$t_L/t_H \left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, tLtH, label=param_text)

        # age
        age = get_age(z, tLtH, age, omega_m, omega_l, H)
        print("age")
        print(age)
        param_text = r"$t/t_H \left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, age, label=param_text)

# 7. intersection probability vs redshift
elif (option == 7):        
    print("Dimensionless Intersection Probability vs Redshift");
    plt.title("Proper Motion Distance vs Redshift");
    plt.ylabel(r"proper motion distance $\Omega_M/\Omega_H$");

    for i in range(len(omegam_array)):
        omega_m = omegam_array[i]
        omega_l = omegal_array[i]
        print(omega_m, omega_l)
        
        dPdz = intersect_prob(z, dPdz, omega_m, omega_l)
        print("dPdz")
        print(dPdz)
        param_text = r"$\left(\Omega_M, \Omega_\Lambda\right) = " + "(" + str(omega_m) + ", " + str(omega_l) + ")$"
        plt.plot(z, dPdz, label=param_text)
    
else:
    print("Error: Please select an integer in the range 1-7.");

plt.axis([zmin, zmax, ymin, ymax]);
plt.legend()
plt.show()
