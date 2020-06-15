# ----------- Import statements ------------

import math;
import numpy as np;
import matplotlib.pyplot as plt;

# ------------ Custom functions ------------

# A-F sums
def A(xlist, ylist, y_uncert):
    A = 0
    for i in range(len(xlist)):
        A += 1.0 / (y_uncert[i])**2
    return A

def B(xlist, ylist, y_uncert):
    B = 0
    for i in range(len(xlist)):
        B += xlist[i] / (y_uncert[i])**2
    return B

def C(xlist, ylist, y_uncert):
    C = 0
    for i in range(len(xlist)):
        C += xlist[i]**2 / (y_uncert[i])**2
    return C

def D(xlist, ylist, y_uncert):
    D = 0
    for i in range(len(xlist)):
        D += xlist[i]**3 / (y_uncert[i])**2
    return D

def E(xlist, ylist, y_uncert):
    E = 0
    for i in range(len(xlist)):
        E += xlist[i]**4 / (y_uncert[i])**2
    return E

def F(xlist, ylist, y_uncert):
    F = 0
    for i in range(len(xlist)):
        F += ylist[i] / (y_uncert[i])**2
    return F
    
def G(xlist, ylist, y_uncert):
    G = 0
    for i in range(len(xlist)):
        G += xlist[i] * ylist[i] / (y_uncert[i])**2
    return G

def H(xlist, ylist, y_uncert):
    H = 0
    for i in range(len(xlist)):
        H += xlist[i]**2 * ylist[i] / (y_uncert[i])**2
    return H

# chi-square
def chi_square(xlist, ylist, y_uncert, b, m, q):
    chi2 = 0
    for i in range(len(xlist)):
        chi2 += (ylist[i] - q*xlist[i]**2 - m*xlist[i] - b)**2 / (y_uncert[i])**2
    return chi2;

# average y-value
def avg(alist):
    avg = 0
    for i in range(len(alist)):
        avg += ylist[i];
    return avg;

# coefficient of determination (r^2)
# r^2 = 1 - SE_ypredicted / SE_yavg (standard errors calculated with respect to y_i)
def r2(xlist, ylist, y_uncert, b, m, q, y_avg):
    r2 = 0
    num = 0
    denom = 0
    se1 = 0
    se2 = 0
    for i in range(len(xlist)):
        se1 = (q*xlist[i]**2 + m*xlist[i] + b - ylist[i])
        num += se1*se1;
        se2 = (ylist[i] - y_avg)
        denom += se2*se2
    r2 = 1.0 - num / denom;
    return r2

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

# Assign the A-F sum values to variables A-H
A = A(xlist, ylist, y_uncert)
B = B(xlist, ylist, y_uncert)
C = C(xlist, ylist, y_uncert)
D = D(xlist, ylist, y_uncert)
E = E(xlist, ylist, y_uncert)
F = F(xlist, ylist, y_uncert)
G = G(xlist, ylist, y_uncert)
H = H(xlist, ylist, y_uncert)

# y = qx^2 + mx + b
denom = -(C*C*C) + (A*C*E) + (2.0*B*C*D) - (A*D*D) - (B*B*E)

b = (
    ( F*(E*C-D*D) + G*(-B*E+C*D) - H*(C*C-B*D) )
     / (denom)
      )
m = (
    ( (F*(C*D-B*E) - G*(C*C-A*E))*(B*C-A*D) + H*(A*A*D*D-2.0*A*B*C*D+B*B*C*C) )
     / ( (B*C - A*D) * denom )
     )
q = (
    ( -F*(C*C-B*D) + G*(B*C-A*D) + H*(A*C-B*B) )
     / (denom)
      )
# Uncertainties in b, m, q
# How do I implement?

# Calculate average y-value
y_avg = avg(ylist);

# Calculate chi-square
chi2 = chi_square(xlist, ylist, y_uncert, b, m, q)
# Calculate degrees of freedom
ndf = len(xlist) - 2        # 2 parameters: a,b
# Calculate closeness of fit (should be as close to 1 as possible)
fit = chi2 / ndf;

# Calculate coefficient of determination r^2 and r
r2 = r2(xlist, ylist, y_uncert, b, m, q, y_avg)
r = math.sqrt(abs(r2))

# ------------ Console output  -------------

# Print the linear regression model
print(r"$y = qx^2 + mx + b$")
print("y = " + str(b) + " + " + str(m) + "x" + str(q) + "x^2")      # equation: y = slope*x + yint
print("b = " + str(b))                          # b value
print("m = " + str(m))                          # m value
print("q = " + str(q))                          # q value
print("chi-square: " + str(chi2))               # chi-square
print("degrees of freedom: " + str(ndf))        # degrees of freedom
print("closeness of fit: " + str(fit))          # closeness of fit (chi-square)
print("coeff. of determ. : r^2 = " + str(r2))   # coefficient of determination
print("correlation: r = " + str(r));            # correlation

# ------------ Plotting output -------------

# Set parameters for plot
print("range of x-values: [" + str(min(xlist)) + " , " + str(max(xlist)) + "]");
print("range of y-values: [" + str(min(ylist)) + " , " + str(max(ylist)) + "]");

xmin = 0.0
xmax = 1.2*max(xlist)
ymin = 0.0
ymax = 1.2*max(ylist)

numdec = 4;                 # number of decimals to round to
# Display the equation and r-value
equation = ( r"$y = " + str(round(q,numdec))
                 + "x^2+" + str(round(m,numdec))
                 + "x+" + str(round(b,numdec)) + "$" 
            + "\n" + r"$r^2=$" + str(r2)
            + "\n" + "ndf = " + str(ndf)
            + "\n" + r"$\chi^2=$" + str(chi2)
            + "\n" + r"$\chi^2/$ndf = " + str(fit)
             )

# Annotate each point
#for x,y in zip(xlist,ylist):
#    label = "(" + "{:.3f}".format(x) + ", " + "{:.3f}".format(y/1e6) + "*10^6)"
#    plt.annotate(label, # this is the text
#                 (x,y), # this is the point to label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,8), # distance from text to points (x,y)
#                 ha='center') # horizontal alignment can be left, right or center


# Plot the axes and labels (need to hardcode xlabel and ylabel)
plt.title(r"$y(x) = qx^2 + mx + b$");
plt.xlabel("x");
plt.ylabel("y");
plt.axis([xmin, xmax, ymin, ymax]);

# Plot the data points
plt.plot(xlist, ylist, 'k.');

plt.errorbar(xlist, ylist, xerr=x_uncert, yerr=y_uncert, fmt='none', ecolor='black');

# Plot the best-fit line
x = np.linspace(xmin,xmax,100);
y = q*x**2 + m*x + b;
plt.plot(x, y, 'black', label=equation);
plt.legend(loc='top right');

plt.show();
