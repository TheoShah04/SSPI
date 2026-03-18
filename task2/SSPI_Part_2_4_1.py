import numpy as np
import matplotlib.pyplot as plt
import csv, sys, os

filepath = "./task2/NASDAQ.csv"

dates, closes = [], []
with open(filepath, "r") as f:
    reader = csv.reader(f)
    header = next(reader)
    for row in reader:
        dates.append(row[0])
        closes.append(float(row[1]))

y_raw = np.array(closes)

y = (y_raw - y_raw.mean()) / y_raw.std()
N = len(y)

# acf
max_order = 10

def acf_biased(y, max_lag):
    N = len(y)
    return np.array([np.sum(y[:N-k] * y[k:]) / N for k in range(max_lag + 1)])

r = acf_biased(y, max_order)

a1_est = r[1] / r[0]

crlb_a1 = (1 - a1_est ** 2) / N
sigma2_est = r[0] * (1 - a1_est**2)

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(a1_range, crlb_vs_a1, 'b-', linewidth=2)
ax.axvline(x=a1_est, color='r', linestyle='--', linewidth=1.5)
ax.axhline(y=crlb_a1, color='r', linestyle=':', alpha=0.5)
ax.plot(a1_est, crlb_a1, 'ro', markersize=8, label=rf'CRLB = {crlb_a1:.2e}, $a_1={a1_est:.4f}$')

ax.set_xlabel(r'$a_1$', fontsize=14)
ax.set_ylabel(r'CRLB', fontsize=14)
ax.set_title(rf'Behaviour of var($\hat{{a}}_1$) as $a_1$ approaches unity', fontsize=14)
ax.legend(fontsize=12)
ax.set_xlim([0, 1])
ax.grid(True, alpha=0.3)