import numpy as np
import matplotlib.pyplot as plt
low = 0
high = 1
N = 1000
x = np.random.uniform(low, high, N)
plt.plot(x, '.', markersize=3)
plt.xlabel("Sample index")
plt.ylabel("Value")
plt.title("Uniformly Distributed Samples")
plt.show()

theoretical_mean = (low + high)/2
sample_mean = x.sum()/N
print("Theoretical Mean: ", theoretical_mean)
print("Sample Mean: ", sample_mean)
print("Mean error: ", theoretical_mean - sample_mean)