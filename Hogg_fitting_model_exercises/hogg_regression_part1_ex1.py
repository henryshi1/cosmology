# ----------- Import statements ------------

import math;
import numpy as np;
import matplotlib.pyplot as plt;

# ------------ Custom functions ------------

# A-F sums
def A(xlist, ylist, y_uncert):
    A = 0;
    for i in range(len(xlist)):
        A += xlist[i] / (y_uncert[i])**2;
    return A;

def B(xlist, ylist, y_uncert):
    B = 0;
    for i in range(len(xlist)):
        B += 1.0 / (y_uncert[i])**2;
    return B;

def C(xlist, ylist, y_uncert):
    C = 0;
    for i in range(len(xlist)):
        C += ylist[i] / (y_uncert[i])**2;
    return C;

def D(xlist, ylist, y_uncert):
    D = 0;
    for i in range(len(xlist)):
        D += (xlist[i])**2 / (y_uncert[i])**2;
    return D;

def E(xlist, ylist, y_uncert):
    E = 0;
    for i in range(len(xlist)):
        E += (xlist[i]) * (ylist[i]) / (y_uncert[i])**2;
    return E;

def F(xlist, ylist, y_uncert):
    F = 0;
    for i in range(len(xlist)):
        F += (ylist[i])**2 / (y_uncert[i])**2;
    return F;

# chi-square
def s_m(xlist, ylist, y_uncert, slope, yint):
    s_m = 0;
    for i in range(len(xlist)):
        s_m += (ylist[i] - slope*xlist[i] - yint)**2 / (y_uncert[i])**2;
    return s_m;

# average y-value
def avg(alist):
    avg = 0;
    for i in range(len(alist)):
        avg += ylist[i];
    return avg;

# coefficient of determination (r^2)
# r^2 = 1 - SE_ypredicted / SE_yavg (standard errors calculated with respect to y_i)
def r2(xlist, ylist, y_uncert, slope, yint, y_avg):
    r2 = 0;
    num = 0;
    denom = 0;
    se1 = 0;
    se2 = 0;
    for i in range(len(xlist)):
        se1 = (slope*xlist[i] + yint - ylist[i])
        num += se1*se1;
        se2 = (ylist[i] - y_avg)
        denom += se2*se2;
    r2 = 1.0 - num / denom;
    return r2;

# ------------------------------------------
#              Enter data here
# ------------------------------------------

# Hardcode these values
xlist = [203,58,210,202,198,158,
         165,201,157,131,166,160,186,125,218,146];
x_uncert = [5,9,4,4,11,
            7,5,5,5,6,6,5,9,8,6,5];

ylist = [495,173,479,504,510,416,
         393,442,317,311,400,337,423,334,533,344];
y_uncert = [21,15,27,14,30,16,
            14,25,52,16,34,31,42,26,16,22];
# -------------- Main program --------------

# Assign the A-F sum values to variables A-F
A = A(xlist, ylist, y_uncert);
B = B(xlist, ylist, y_uncert);
C = C(xlist, ylist, y_uncert);
D = D(xlist, ylist, y_uncert);
E = E(xlist, ylist, y_uncert);
F = F(xlist, ylist, y_uncert);

# y = ax + b
slope = (B*E - A*C) / (B*D - A*A)
yint = (C*D - A*E) / (B*D - A*A)
# Uncertainties in slope and y-intercept
sigma_slope2 = B / (B*D-A*A)
sigma_slope = math.sqrt(sigma_slope2)
sigma_yint2 = D / (B*D-A*A)
sigma_yint = math.sqrt(sigma_yint2)

# Calculate average y-value
y_avg = avg(ylist);

# Calculate chi-square
s_m = s_m(xlist, ylist, y_uncert, slope, yint);
# Calculate degrees of freedom
ndf = len(xlist) - 2;        # 2 parameters: a,b
# Calculate closeness of fit (should be as close to 1 as possible)
fit = s_m / ndf;

# Calculate coefficient of determination r^2 and r
r2 = r2(xlist, ylist, y_uncert, slope, yint, y_avg);
if (slope > 0):
    r = math.sqrt(r2);
else:
    r = -math.sqrt(r2);

# ------------ Console output  -------------

# Print the linear regression model
print("y = " + str(yint) + " + " + str(slope) + "x" );      # equation: y = slope*x + yint
print("y-intercept = " + str(yint));                        # yint value
print("uncertainty in y-intercept = " + str(sigma_yint));   # yint uncertainty
print("slope = " + str(slope));                             # slope value
print("variance in slope = " + str(sigma_slope2));          # slope variance
print("uncertainty in slope = " + str(sigma_slope));        # slope uncertainty
print("chi-square: " + str(s_m));                           # chi-square
print("degrees of freedom: " + str(ndf));                   # degrees of freedom
print("closeness of fit: S_m/ndf = " + str(fit));           # closeness of fit (chi-square)
print("coeff. of determ. : r^2 = " + str(r2));              # coefficient of determination
print("r = " + str(r));                                     # correlation

# ------------ Plotting output -------------

# Set parameters for plot
print("range of x-values: [" + str(min(xlist)) + " , " + str(max(xlist)) + "]");
print("range of y-values: [" + str(min(ylist)) + " , " + str(max(ylist)) + "]");

print("Enter min and max for axes of plot:");
xmin = 0.0#float(input("xmin: "));
xmax = 1.2*max(xlist)#float(input("xmax: "));
ymin = 0.0#float(input("ymin: "));
ymax = 1.2*max(ylist)#float(input("ymax: "));

numdecimals = 5;
# Display the equation and r-value
equation = ( "y = " + str(slope) + "x" + " + " + str(yint)
            + "\n" + "\u03C3_slope = " + str(sigma_slope)
            + "\n" + "\u03C3_slope^2 = " + str(sigma_slope2)
            + "\n" + "r^2 = " + str(r2)
            + "\n" + "ndf = " + str(ndf)
            + "\n" + "\u03C7^2 = " + str(s_m)
            #+ "\n" + "\u03C7^2/ndf = " + str(fit)
            + "\n" + "S_m/ndf = " + str(fit) );

# Annotate each point
#for x,y in zip(xlist,ylist):
#    label = "(" + "{:.3f}".format(x) + ", " + "{:.3f}".format(y/1e6) + "*10^6)"
#    plt.annotate(label, # this is the text
#                 (x,y), # this is the point to label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,8), # distance from text to points (x,y)
#                 ha='center') # horizontal alignment can be left, right or center


# Plot the axes and labels (need to hardcode xlabel and ylabel)
plt.title("y(x)");
plt.xlabel("x");
plt.ylabel("y");
plt.axis([xmin, xmax, ymin, ymax]);

# Plot the data points
plt.plot(xlist, ylist, 'k.');

plt.errorbar(xlist, ylist, xerr=x_uncert, yerr=y_uncert, fmt='none', ecolor='black');

# Plot the best-fit line
x = np.linspace(xmin,xmax,100);
y = slope*x + yint;
plt.plot(x, y, 'black', label=equation);
plt.legend(loc='top right');

plt.show();
