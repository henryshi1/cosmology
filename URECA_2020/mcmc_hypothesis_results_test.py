import numpy as np
import scipy.stats as st

# alpha level
alpha_low = 0.025
alpha_high = 1.0-0.025

# z-score for omegam-values
omegam_exp = -0.337   # experimental omega_m
sigma_omegam_exp = 0.066
omegam_thr = 0.307   # theoretical omega_m
sigma_omegam_thr = 0.012

z_omegam = (omegam_exp - omegam_thr) / np.sqrt(sigma_omegam_exp**2 + sigma_omegam_thr**2)
p_omegam = st.norm.cdf(z_omegam)

# z-score for w-values
wexp = -1.182   # experimental w
sigma_wexp = 0.521
wthr = -1.026   # theoretical w
sigma_wthr = 0.041

z_w = (wexp - wthr) / np.sqrt(sigma_wexp**2 + sigma_wthr**2)
p_w = st.norm.cdf(z_w)

print("Assuming a two-tailed alpha=0.05 significance level, which means that a given p-value must be less than 0.025 or greater than 1-0.025=0.975 in order for result to be significant.")
print("\n-------------------------------\n")

print("z-score for omega_m: " + str(z_omegam))
print("p-value for omega_m: " + str(p_omegam))
if (p_omegam < alpha_low or p_omegam > alpha_high):
    print("omega_m is significantly different than hypothesized")
else:
    print("omega_m is not significantly different than hypothesized")
print("- - - - - - - - - - - - - - - -")
print("z-score for w: " + str(z_w))
print("p-value for w: " + str(p_w))
if (p_w < alpha_low or p_w > alpha_high):
    print("w is significantly different than hypothesized")
else:
    print("w is not significantly different than hypothesized")
