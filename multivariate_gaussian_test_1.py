import scipy.stats as stats
import numpy as np

# Initialize the x- and y-arrays that we will use to plot
x_axis = []
mean = []
variance = []
skew_list = []
kurt_list = []

# generate gaussian pdf of x
#mu = 2.0
#variance = 2.0
#sigma = math.sqrt(variance)
#x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
#plot(x, stats.norm.pdf(x, mu, sigma))
#show()

# proposal function q(x'|x)
def q(x):
    # draw from a multivariate gaussian with mean [0,0] and covariance matrix I
    # Add the randomly drawn number to x and return the sum
    mu = np.array([0.0, 0.0])
    sigma = np.matrix([ [1.0,0.0], [0.0,1.0] ])
    q = stats.multivariate_normal.rvs(mean =mu, cov=sigma)
    return x + q

x = np.array([0.0, 0.0])
x1 = q(x)
print(q(x))
