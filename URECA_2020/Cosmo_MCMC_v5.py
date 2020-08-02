# --------------------------------- Import everything ---------------------------------------

import scipy.stats as stats
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import csv
import time

# ------------------------------- Collect data from file ------------------------------------

# Begin by collecting csv data and making a scatter plot of it
zs = []
dm_obs = []
sigmas = []

with open('ps1_data.txt', newline='') as csvfile:
    ps1_data = csv.reader(csvfile, delimiter=' ')
    i = 0
    for row in ps1_data:
        if (i > 0):
            zs.append(float(row[1]))
            dm_obs.append(float(row[4]))
            sigmas.append(float(row[5]))
        i += 1
csvfile.close()

# -------------------------------- Variable Dictionary --------------------------------------

# Constants
omega_k = 0  # critical density
c = 2.9979e5  # speed of light (km/s)
omega_l = None  # energy density
pressure = 0.0  # pressure of the universe

# Parameters we are trying to find
omega_m = None  # mass density parameter
w = None  # equation of state
dm_offset = None  # offset of the model from the data

# Variable dictionary
H = None  # Hubble constant (km/s/Mpc)
z = None  # redshift
D_C = None  # comoving distance, line-of-sight (m)
D_M = None  # comoving distance, transverse (m)
D_A = None  # angular diameter distance (m)
D_L = None  # luminosity distance (m)
V_C = None  # comoving volume (m^3)
D_H = None  # Hubble Distance (Mpc)
array_size = None  # size of each array, which is the number of z values
param_text = ""  # label for the plots
npts = None  # Number of points (omega_m, w) we use

# Calculate Hubble distance
H = 70.0  # Hubble constant
D_H = c / H  # Hubble distance

# ----------------------------------- Function Dictionary -----------------------------------

# --------------------------------- E(z) and DCDH integral ----------------------------------

# Calculate E(z), which in this situation depends on w
def E(z_i, w, omega_m, omega_l):
    omega_k = 1.0 - omega_m - omega_l
    # Assume a flat universe, i.e. omega_k = 0
    omega_k = 0.0

    inside = omega_m * (1.0 + z_i) ** 3 + omega_k * (1.0 + z_i) ** 2 + omega_l * (1.0 + z_i) ** (3.0 * (1.0 + w))

    if (inside < 0):
        print("inside = " + str(inside))
    #inside = abs(inside)
    E = np.sqrt(inside)
    return E

# Calculate the integral D_C/D_H, store it in an array
def DCDH_int(z_limit, w, omega_m, omega_l):
    f = lambda z_i: 1.0 / E(z_i, w, omega_m, omega_l)
    i = integrate.quad(f, 0, z_limit)
    return i[0]

# ----------------------------- Proper Motion distance integrals -----------------------------

# proper motion distance for omega_l == 0
def prop_motion_0(z, w, omega_m, omega_l):
    # find omega_k = 1 - omega_m - omega_l
    # omega_k = 1.0 - omega_m - omega_l
    # Continuing our flat universe assumption, we have
    omega_k = 0.0

    # Calculate D_M/D_H for different universe geometries
    DCDH = DCDH_int(z, w, omega_m, omega_l)

    # open universe (omega_k > 0)
    if (omega_k > 1.0e-6):
        return 1 / np.sqrt(omega_k) * np.sinh(np.sqrt(omega_k) * DCDH)
    # closed universe (omega_k < 0)
    elif (omega_k < -1.0e-6):
        return 1 / np.sqrt(abs(omega_k)) * np.sin(np.sqrt(abs(omega_k)) * DCDH)
    # flat universe (omega_k = 0)
    else:
        return 2.0 * (2.0 - omega_m * (1.0 - z) - (2 - omega_m) * np.sqrt(1.0 + omega_m * z)) / (
                    omega_m * omega_m * (1.0 + z))

# proper motion distance for omega_l != 0
def prop_motion(z, w, omega_m, omega_l):
    # find omega_k = 1 - omega_m - omega_l
    omega_k = 1.0 - omega_m - omega_l
    # Continuing our flat universe assumption, we have
    omega_k = 0.0

    # Calculate D_C/D_H integral
    DCDH = DCDH_int(z, w, omega_m, omega_l)

    # Calculate D_M/D_H for different universe geometries
    # open universe
    if (omega_k > 1.0e-6):
        return 1 / np.sqrt(omega_k) * np.sinh(np.sqrt(omega_k) * DCDH)
    # closed universe
    elif (omega_k < -1.0e-6):
        return 1 / np.sqrt(abs(omega_k)) * np.sin(np.sqrt(abs(omega_k)) * DCDH)
    # flat universe
    else:
        return DCDH

# --------------------------------------------------------------------------------------------

# distance modulus
def dist_mod(z, w, D_H, omega_m, omega_l, DMDH):
    # initialize the distance modulus value
    mu = 0.0

    # convert D_H from Mpc to pc
    DH = D_H * 1.0e6

    # Determine D_M/D_H
    if (omega_l == 0):
        DMDH = prop_motion_0(z, w, omega_m, omega_l)
    else:
        DMDH = prop_motion(z, w, omega_m, omega_l)

    # calculate each value of mu array
    # for i in range(len(z)):
    # mu[i] = 5.0 * ( np.log10(1.0+z[i]) + np.log10(DCDH_int(z[i], omega_m, omega_l)) + np.log10(D_H/10) )
    mu = 5.0 * (np.log10(1.0 + z) + np.log10(DMDH) + np.log10(D_H / 10))

    return mu

# --------------------------------------------------------------------------------------------

# ----------------------- Proposal and log likelihood functions defined ----------------------

# ----------------------- log likelihood function -----------------------

# density function lnf(x) is the log likelihood of the data given parameters
# takes a vector x = (x1, x2, x3) as input
def make_log_likelihood(zs, dm_obs, sigmas, D_H):
    def log_likelihood(omega_m, w, dm_offset):
        # Assume a flat universe, i.e. omega_l = 1 - omega_m, or omega_k = 0
        omega_l = 1.0 - omega_m

        # initialize the log likelihood function
        ln_likelihood = 0.0

        # loop for all values of z in the z-array
        for i in range(len(zs)):
            # initialize for each iteration
            z_i = zs[i]
            sigma_dm_i = sigmas[i]
            dm_obs_i = dm_obs[i]        # observed distance modulus for the given z-value

            # calculate predicted value of distance modulus based on z-value and parameters
            # predicted = [distance modulus calculation for given z] + offset
            DMDH_i = prop_motion(z_i, w, omega_m, omega_l)  # proper motion distance
            dm_pred_i = dist_mod(z_i, w, D_H, omega_m, omega_l, DMDH_i) + dm_offset  # distance modulus

            # calculate each term
            term = np.log(1.0 / (sigma_dm_i * np.sqrt(2.0 * np.pi))) \
                   - 0.5 * (dm_obs_i - dm_pred_i) * (dm_obs_i - dm_pred_i) / (sigma_dm_i * sigma_dm_i)

            # increment log-likelihood
            ln_likelihood += term

            # print statement to debug code
            # print("x: " + str(xi) + " | sigma: " + str(sigma_i) + " | Term being added: " + str(term))

        return ln_likelihood

    return log_likelihood

# -------------------------------- proposal function --------------------------------

# proposal function q(x'|x) where x is a 3D vector. Proposal has a mean of x and each component has variance 1
def q(x, num_iter, zs, dm_obs, sigmas, D_H):
    # --------------------------------  Variables --------------------------------
    num_iter = num_iter  # num_iter = number of iterations
    vec_size = len(x)  # vec_size = size of each vector (number of dimensions)
    xchain = [None] * num_iter  # Initialize the Markov chain of x-vectors

    # propose num_iter matter density values (omega_m) in a uniform range 0 to 1
    # resulting in a chain of omega_m values of length num_iter
    #omega_ms = np.random.uniform(0.0, 1.0, num_iter)
    omega_ms = np.empty(num_iter)
    omega_ms.fill(0.337)

    w = 0.0  # w, Initialize as 0
    dm_offset = 0.0  # distance modulus offset, Initialize as 0

    # temporary vectors for x
    old = x  # old x-vector, initialized as x
    new = None  # new x-vector

    # log-likelihood function
    lnf = make_log_likelihood(zs, dm_obs, sigmas, D_H)

    # -------------------------------- loop --------------------------------

    for i in range(num_iter):  # q is a 3D vector, so 3 components

        # initialize the "old" vector
        if (i > 0):
            old = xchain[i - 1]

        # ----- Print statement for debugging -----
        print("Iteration number: i = " + str(i))
        print("old x: " + str(old))
        # -----------------------------------------

        omega_m = 0.337
        w = -1.182
        dm_offset = old[2]

        omega_m_new = 0.337
        w_new = -1.182

        # propose dm_offset from a normal distribution with mean dm_offset (the old value) and variance 1
        mu = dm_offset
        sigma = 1.0
        dm_offset_new = np.random.normal(mu, sigma, 1)[0]

        # alternatively propose dm_offset_new
        #dm_offset_new = np.random.uniform(-30,30)

        # generate each new entry of xchain. The components of each entry are as follows:
        new = [omega_m_new, w_new, dm_offset_new]

        # ----- Print statement for debugging -----
        print("new x: " + str(new))
        # -----------------------------------------

        # Conditional statement: if lnfnew - lnfold > benchmark then tack on new value
        # otherwise, tack on old array copy the old value into the new value
        lnfnew = lnf(omega_m_new, w_new, dm_offset_new)
        lnfold = lnf(omega_m, w, dm_offset)
        bench_in = np.random.uniform(0.0,1.0)
        benchmark = np.log(bench_in)
        #benchmark = np.random.uniform(0.0,1.0)

        # ----- Print statement for debugging -----
        #print("lnfnew: " + str(lnfnew))
        #print("lnfold: " + str(lnfold))
        #print("lnfnew - lnfold = " + str(lnfnew - lnfold))
        #print("benchmark inside: " + str(bench_in))
        #print("benchmark: " + str(benchmark))
        # -----------------------------------------

        if (lnfnew - lnfold > benchmark):   # benchmark inside term is from 0 to 1,
            xchain[i] = new                 # so benchmark itself is from -inf to 0.
            print("new chain appended")     # lnf() flips the sign of inside!
        else:                               # Therefore, the inequality sign should be flipped
            xchain[i] = old
            print("old chain appended")

        # dummy statement, consider removing
        #xchain[i] = new

        print()
        print("------------------------------------------")
        #time.sleep(0.0)

    return xchain

# --------------------------------------- Driver function ---------------------------------------
# test q

x = [0.337, -1.182, 0.0]  # seed vector
print("x before proposal: " + str(x))

# Generate Markov chain of x-vectors
num_iterations = 10000                     # number of iterations
xchain = q(x, num_iterations, zs, dm_obs, sigmas, D_H)

# The final x-vector is the one we want
xfinal = xchain[-1]
print("xinitial = " + str(x))
print("xfinal = " + str(xfinal))
print("xchain: " + str(xchain))
print("length of xchain: " + str(len(xchain)))

# --------------------------------------- Display results ---------------------------------------
# Extract V- and R-coordinates from the array of vectors
omegams = []    # list of omega_m values
ws = []         # list of w values
offsets = []    # list of dm offsets

j = 0
burnin_point = -1
for i in xchain:
    if (j > burnin_point):
        omegams.append(i[0])
        ws.append(i[1])
        offsets.append(i[2])
    j += 1
    #print("length of each array: " + str(len(omegams)))
    # Variation of the loop: use an if statement to omit burn-in points

# Print the average of each parameter
print()
print("------------------------------------------")
print("average omega_m = " + str(np.mean(omegams)))
print("uncertainty in omega_m = "+ str(np.std(omegams)))
print("average w = " + str(np.mean(ws)))
print("uncertainty in w = "+ str(np.std(ws)))
print("average dm_offset = " + str(np.mean(offsets)))
print("uncertainty in dm_offset = "+ str(np.std(offsets)))

# Trace the evolution of omega_m, w, dm_offset values
plt.subplot(321)
plt.title('omega_m values')
plt.plot(omegams)

plt.subplot(323)
plt.title('w values')
plt.plot(ws)

plt.subplot(325)
plt.title('offset values')
plt.plot(offsets)

# Plot histograms of each  parameter
plt.subplot(322)
plt.hist(omegams)

plt.subplot(324)
plt.hist(ws)

plt.subplot(326)
plt.hist(offsets)

plt.tight_layout()

plt.show()
