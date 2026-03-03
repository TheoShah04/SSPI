import numpy as np
import matplotlib.pyplot as plt
low = 0
high = 1
N = 1000
x = np.random.uniform(low, high, N)
theoretical_mean = (low + high)/2
sample_mean = x.sum()/N

theoretical_SD = np.sqrt(np.mean(((x - 0.5) ** 2)))
sample_SD = np.sqrt(1/(N-1) * (np.sum(np.square(x-sample_mean))))
print("Theoretical SD: ", theoretical_SD)
print("Sample SD: ", sample_SD)
print("SD error: ", sample_SD - theoretical_SD)
