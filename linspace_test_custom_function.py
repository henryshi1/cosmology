import numpy as np
import math

def prop_motion_0(z, DMDH, density_m, density_l):
    DMDH = np.copy(z)
    DMDH = 2.0*( 2.0 - density_m*(1.0-z) - (2-density_m)*np.sqrt(1.0+density_m*z) ) / ( density_m*density_m*(1.0+z) );
    print(len(DMDH))
    print(DMDH)
    return

z = np.linspace(0,5,100)
DMDH = np.array([])
density_m = 0.3
density_l = 0.7

prop_motion_0(z, DMDH, density_m, density_l)
print(DMDH)
