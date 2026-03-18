import numpy as np
import matplotlib.pyplot as plt

def rp3(M, N):
    a = 0.5 
    m = 3.0

    v = (np.random.rand(M, N) - 0.5) * m + a # All realisations are random. Same distribution across time and realisations. Ergodic process.
    return v

def pdf(v, bins=50):
    counts, edges = np.histogram(v, bins=bins, density=False)
    dx = edges[1] - edges[0]
    p = counts / (np.sum(counts) * dx)
    x = edges[:-1]
    return p, x

N = 1000
M = 100
N_list = [100, 1000, 10000]
process_3 = rp3(M, N)

a = 0.5
m = 3.0

x_th = np.linspace(a - m/2, a + m/2, 500)
pdf_th = np.ones_like(x_th) / m

for N in N_list:
    process_3 = rp3(M, N)
    p, x = pdf(process_3, bins=50)

    plt.figure()
    plt.bar(x, p, width=x[1] - x[0], alpha=0.6,
            label=f"Estimated PDF (N={N})")
    plt.plot(x_th, pdf_th, 'r', linewidth=2,
             label="Theoretical PDF")

    plt.xlabel("x")
    plt.ylabel("Probability density")
    plt.title(f"PDF estimate for rp3 with N = {N}")
    plt.legend()
    plt.grid(True)
