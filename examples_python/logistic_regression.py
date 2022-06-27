import matplotlib.pyplot as plt
import numpy as np


# design parameters
beta0 = - 20.0
beta1 = 0.7
beta1 = 0.9
fmax = 40.0

beta0 = - 5.0
beta1 = 10.0
fmax = 1.0

def logistic_regresson(beta0, beta1, x):
    return 1.0 / (1.0 + np.exp(-beta1 * x - beta0))

x = np.linspace(0.0, fmax, 100)
y = logistic_regresson(beta0, beta1, x)

plt.plot(x, y)
plt.show()