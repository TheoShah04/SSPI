import numpy as np
import matplotlib.pyplot as plt
low = 0
high = 1
M = 10
N = 1000
x = np.random.uniform(low, high, size=(N, M))
theoretical_mean = (low + high)/2
sample_mean = np.sum(x, axis=0)/N
# theoretical_SD = (high-low)/np.sqrt(12)
theoretical_SD = np.sqrt(np.mean(((x - 0.5) ** 2)))
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

