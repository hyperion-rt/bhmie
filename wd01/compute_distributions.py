import numpy as np
from numpy import pi
from scipy.special import erf

A0 = [3.5e-4, 30.e-4]  # microns
RHO = 2.24  # cm^2/g
BCF = [0.75, 0.25]
SIGMA = 0.4
MC = 1.99442341e-23  # g


def F(a, beta, a_t):
    if beta >= 0:
        return 1. + beta * a / a_t
    else:
        return 1. / (1 - beta * a / a_t)


def D(a, b_c):

    B = []
    for i in range(2):
        term_1 = 3. / (2 * pi) ** 1.5
        term_2 = np.exp(-4.5 * SIGMA ** 2) \
               / (RHO * (A0[i] * 1.e-4) ** 3 * SIGMA)
        term_3 = (BCF[i] * b_c * MC) \
               / (1. + erf(3 * SIGMA / np.sqrt(2) \
                           + np.log(A0[i] / 3.5e-4) / (SIGMA * np.sqrt(2))))
        B.append(term_1 * term_2 * term_3)

    array = np.zeros(a.shape)
    for i in range(2):
        array += B[i] / a * np.exp(- 0.5 * (np.log(a / A0[i]) / SIGMA) ** 2.)
    return array


def dn_over_da(a, C, a_t, a_c, alpha, beta, b_c):
    n = np.zeros(a.shape)
    n[(a > 3.5e-4) & (a < a_t)] = 1.
    n[a >= a_t] = np.exp(-((a[a >= a_t] - a_t) / a_c) ** 3)
    n *= C / a * (a / a_t) ** alpha * F(a, beta, a_t)
    n[a >= 3.5e-4] += D(a[a >= 3.5e-4], b_c)
    return n

# Generate grain sizes
n_a = 1000
amin = 1.e-4  # microns
amax = 100.  # microns
a = np.logspace(np.log10(amin), np.log10(amax), n_a)

# Open parameters for different dust types
f = open('data/table_1', 'rb')

# Loop through them at for each, create the size distributions
for line in f.readlines():

    cols = line.split()

    # Main parameters
    r_v = float(cols[0])
    b_c = float(cols[1])
    case = cols[2]

    # Graphite
    alpha, beta, a_t, a_c, C = [float(col) for col in cols[3:8]]
    n = dn_over_da(a, C, a_t, a_c, alpha, beta, b_c * 1.e-5)
    np.savetxt('sizes/wd01_graphite_%3.1f_%3.1f_%1s' \
               % (r_v, b_c, case), zip(a, n), fmt="%12.4e")

    # Silicate
    alpha, beta, a_t, C = [float(col) for col in cols[8:12]]
    a_c = 0.1
    n = dn_over_da(a, C, a_t, a_c, alpha, beta, 0.)
    np.savetxt('sizes/wd01_silicate_%3.1f_%3.1f_%1s' \
               % (r_v, b_c, case), zip(a, n), fmt="%12.4e")

f.close()
