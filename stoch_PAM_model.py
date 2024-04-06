import numpy as np
import matplotlib.pyplot as plt

# Generate noise again
ns1 = np.random.uniform(-1, 1, 500000)
ns2 = np.random.uniform(-1, 1, 500000)

# Allocate memory for variables again
m_til1 = np.zeros_like(ns1)
m_til2 = np.zeros_like(ns2)
u1 = np.zeros_like(ns1)
u2 = np.zeros_like(ns2)
m1 = np.zeros_like(ns1)
m2 = np.zeros_like(ns2)

# Initial conditions
m_til1[0], m_til1[1] = -1., 0.5
m_til2[0], m_til2[1] = 1., 0.5
u1[0], u2[0] = 10., 10.

# Constants
F1, F2 = 1./9., 1./8.
b11, b12, b21, b22 = 0.04, 0.052, -0.025, 0.00

# Update m_tilda, u, and m iteratively again
for t in range(2, len(ns1)):
    m_til1[t] = 0.6 * m_til1[t - 1] - 0.3 * m_til1[t - 2] + ns1[t]
    m_til2[t] = 0.6 * m_til2[t - 1] - 0.3 * m_til2[t - 2] + ns2[t]

for t in range(1, len(ns1)):
    m1[t] = m_til1[t] + b11 * u1[t - 1] + b12 * u2[t - 1]
    m2[t] = m_til2[t] + b21 * u1[t - 1] + b22 * u2[t - 1]
    u1[t] = u1[t - 1] + (m1[t] - F1 * u1[t - 1])
    u2[t] = u2[t - 1] + (m2[t] - F2 * u2[t - 1])


# Cross correlation using numpy's correlate function for simplicity
mxlag = 60
lags = np.arange(-mxlag, mxlag+1, 1)

def cross_correlation(x, y, maxlag):
    """Simple cross correlation function."""
    corr = np.correlate(x, y, mode='full')
    start = len(x) - 1 - maxlag
    end = len(x) + maxlag
    return corr[start:end]
    
# Calculating cross-correlations again
ccr2 = cross_correlation(u1, u1, mxlag)
ccr3 = cross_correlation(u2, u2, mxlag)
ccr4 = cross_correlation(u1, u2, mxlag)
ccr5 = cross_correlation(m1, m2, mxlag)
ccr6 = cross_correlation(u1, m1, mxlag)
ccr7 = cross_correlation(u2, m2, mxlag)
ccr8 = cross_correlation(u2, m1, mxlag)
ccr9 = cross_correlation(u1, m2, mxlag)

# Set up the plots again
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

axs[0, 0].plot(lags, ccr2, label='u1u1')
axs[0, 1].plot(lags, ccr3, label='u2u2')
axs[0, 2].plot(lags, ccr4, label='u1u2')
axs[1, 0].plot(lags, ccr5, label='m1m2')
axs[1, 1].plot(lags, ccr6, label='u1m1')
axs[1, 2].plot(lags, ccr7, label='u2m2')

for ax in axs.flat:
    ax.legend()
    ax.grid(True)
    ax.set(xlabel='Lag', ylabel='Cross-correlation')

plt.tight_layout()

# Save the plot to a file
plt.savefig("cross_correlation_plots.png")
