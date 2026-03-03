import numpy as np
import matplotlib.pyplot as plt
N = 1000
x = np.random.randn(N)
plt.plot(x, '.', markersize=3)
plt.xlabel("Sample index")
plt.ylabel("Value")
plt.title("Uniformly Distributed Samples")
plt.show()

theoretical_mean = 0
sample_mean = x.sum()/N
print("Theoretical Mean: ", theoretical_mean)
print("Sample Mean: ", sample_mean)
print("Mean error: ", theoretical_mean - sample_mean)

theoretical_SD = 1
sample_SD = np.sqrt(1/(N-1) * (np.sum(np.square(x-sample_mean))))
print("Theoretical SD: ", theoretical_SD)
print("Sample SD: ", sample_SD)
print("SD error: ", sample_SD - theoretical_SD)

import random as rng
M = 10
N = 1000
rng = np.random.default_rng()
x = rng.standard_normal((N, M))
theoretical_mean = 0
sample_mean = np.sum(x, axis=0)/N
theoretical_SD = 1
sample_SD = np.sqrt(1/(N-1) * (np.sum(np.square(x-sample_mean), axis=0)))
plt.plot(sample_mean, '.', markersize=3)
plt.axhline(theoretical_mean, color='red', linestyle='--', linewidth=1)
plt.xlabel("Realisation Index")
plt.ylabel("Sample mean")
plt.title("Sample means over 10 realisations")
plt.show()

plt.plot(sample_SD, '.', markersize=3)
plt.axhline(theoretical_SD, color='red', linestyle='--', linewidth=1)
plt.xlabel("Realisation Index")
plt.ylabel("Sample standard deviation")
plt.title("Sample standard deviation over 10 realisations")
plt.show()

N = 100000
x_hist = np.random.randn(N)
x = np.linspace(-4, 4, 1000)
theoretical_pdf = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x**2)
plt.hist(x_hist, bins=50, density=True, color='blue', edgecolor='black')
plt.plot(x, theoretical_pdf, 'r--', linewidth=2)
plt.xlabel('Value')
plt.ylabel('pdf')
plt.title('Theoretical and Estimate pdf')
plt.show()