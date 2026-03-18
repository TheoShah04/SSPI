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

M = 100
N = 100

def generateSamplesAndPlots(M, N):
    v_1 = rp1(M, N)
    v_2 = rp2(M, N)
    v_3 = rp3(M, N)
    ensemble_mean_v_1 = np.mean(v_1, axis=0)
    ensemble_SD_v_1 = np.sqrt(1/(M-1) * (np.sum(np.square(v_1-ensemble_mean_v_1), axis=0)))
    ensemble_mean_v_2 = np.mean(v_2, axis=0)
    ensemble_SD_v_2 = np.sqrt(1/(M-1) * (np.sum(np.square(v_2-ensemble_mean_v_2), axis=0)))
    ensemble_mean_v_3 = np.mean(v_3, axis=0)
    ensemble_SD_v_3 = np.sqrt(1/(M-1) * (np.sum(np.square(v_3-ensemble_mean_v_3), axis=0)))

    plt.plot(ensemble_mean_v_1, '.', markersize=3)
    plt.xlabel("Time Index")
    plt.ylabel("Ensemble mean")
    plt.title("Ensemble means over 100 time indexes (Process 1)")
    plt.show()

    plt.plot(ensemble_SD_v_1, '.', markersize=3)
    plt.xlabel("Time Index")
    plt.ylabel("Ensemble SD")
    plt.title("Ensemble SD over 100 time indexes (Process 1)")
    plt.show()

    plt.plot(ensemble_mean_v_2, '.', markersize=3)
    plt.xlabel("Time Index")
    plt.ylabel("Ensemble mean")
    plt.title("Ensemble means over 100 time indexes (Process 2)")
    plt.show()

    plt.plot(ensemble_SD_v_2, '.', markersize=3)
    plt.xlabel("Time Index")
    plt.ylabel("Ensemble SD")
    plt.title("Ensemble SD over 100 time indexes (Process 2)")
    plt.show()

    plt.plot(ensemble_mean_v_3, '.', markersize=3)
    plt.xlabel("Time Index")
    plt.ylabel("Ensemble mean")
    plt.title("Ensemble means over 100 time indexes (Process 3)")
    plt.show()

    plt.plot(ensemble_SD_v_3, '.', markersize=3)
    plt.xlabel("Time Index")
    plt.ylabel("Ensemble SD")
    plt.title("Ensemble SD over 100 time indexes (Process 3)")
    plt.show()

generateSamplesAndPlots(M, N)