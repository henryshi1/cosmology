import numpy as np
import scipy.stats as stats
x = np.linspace(0, 5, 10, endpoint=False)
print(x)
y = stats.multivariate_normal.pdf(x, mean=2.5, cov=0.5)
print(y)
