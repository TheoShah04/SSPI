import numpy as np
import matplotlib.pyplot as plt
low = 0
high = 1
M = 10
N = 1000
x_hist = np.random.uniform(low, high, N*100)
theoretical_pdf = 1
plt.hist(x_hist, bins=15, density=True, color='blue', edgecolor='black')
plt.axhline(theoretical_pdf, color='red', linestyle='--', linewidth=1)
plt.xlabel('Value')
plt.ylabel('pdf')
plt.title('Theoretical and Estimate pdf')
plt.show()
