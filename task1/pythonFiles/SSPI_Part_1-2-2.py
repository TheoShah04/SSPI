import numpy as np
import matplotlib.pyplot as plt

def rp1(M, N):
    a = 0.02
    b = 5.0

    Mc = np.ones((M, 1)) * b * np.sin(np.arange(1, N + 1) * np.pi / N)
    Ac = a * np.ones((M, 1)) * np.arange(1, N + 1) # Offset

    v = (np.random.rand(M, N) - 0.5) * Mc + Ac # Zero mean noise applied to deterministic sin wave
    return v

def rp2(M, N): # All rows have the same value
    Ar = np.random.rand(M, 1) * np.ones((1, N)) # Constant in offset for each realisation.
    Mr = np.random.rand(M, 1) * np.ones((1, N)) # Constant in amplitude for each realisation.

    v = (np.random.rand(M, N) - 0.5) * Mr + Ar # Zero mean noise applied to non-deterministic noise
    return v

def rp3(M, N):
    a = 0.5 
    m = 3.0

    v = (np.random.rand(M, N) - 0.5) * m + a # All realisations are random. Same distribution across time and realisations. Ergodic process.
    return v

def generateSamplesAndPlots_timeAverage(M, N):
    v_1 = rp1(M, N)
    v_2 = rp2(M, N)
    v_3 = rp3(M, N)

    time_mean_v1 = np.mean(v_1, axis=1)
    time_std_v1  = np.std(v_1, axis=1, ddof=1)

    time_mean_v2 = np.mean(v_2, axis=1)
    time_std_v2  = np.std(v_2, axis=1, ddof=1)

    time_mean_v3 = np.mean(v_3, axis=1)
    time_std_v3  = np.std(v_3, axis=1, ddof=1)

    plt.plot(time_mean_v1, '.', markersize=6)
    plt.xlabel("Realization Index")
    plt.ylabel("Time Mean")
    plt.title("Time Mean of Each Realization (Process 1)")
    plt.show()

    plt.plot(time_std_v1, '.', markersize=6)
    plt.xlabel("Realization Index")
    plt.ylabel("Time SD")
    plt.title("Time Standard Deviation of Each Realization (Process 1)")
    plt.show()

    plt.plot(time_mean_v2, '.', markersize=6)
    plt.xlabel("Realization Index")
    plt.ylabel("Time Mean")
    plt.title("Time Mean of Each Realization (Process 2)")
    plt.show()

    plt.plot(time_std_v2, '.', markersize=6)
    plt.xlabel("Realization Index")
    plt.ylabel("Time SD")
    plt.title("Time Standard Deviation of Each Realization (Process 2)")
    plt.show()

    plt.plot(time_mean_v3, '.', markersize=6)
    plt.xlabel("Realization Index")
    plt.ylabel("Time Mean")
    plt.title("Time Mean of Each Realization (Process 3)")
    plt.show()

    plt.plot(time_std_v3, '.', markersize=6)
    plt.xlabel("Realization Index")
    plt.ylabel("Time SD")
    plt.title("Time Standard Deviation of Each Realization (Process 3)")
    plt.show()

M = 4
N = 1000
generateSamplesAndPlots_timeAverage(M, N)
