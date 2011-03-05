import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Opacity to extinction

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for model in ['3.1_6.0', '4.0_4.0', '5.5_3.0']:

    orig = np.loadtxt('d03_%s_A.orig' % model)
    ax.loglog(orig[:, 0], orig[:, 3], color='black', lw=1)

    new = np.loadtxt('d03_%s_A.summary' % model)
    ax.loglog(new[:, 0], new[:, 3] * 1.67e-24, color='red', lw=1)

fig.savefig('comparison_chi.png')

# Average scattering angle

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for model in ['3.1_6.0', '4.0_4.0', '5.5_3.0']:

    orig = np.loadtxt('d03_%s_A.orig' % model)
    ax.plot(orig[:, 0], orig[:, 2], color='black', lw=1)

    new = np.loadtxt('d03_%s_A.summary' % model)
    ax.plot(new[:, 0], new[:, 4], color='red', lw=1)

ax.set_xscale('log')
fig.savefig('comparison_g.png')

# Albedo

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for model in ['3.1_6.0', '4.0_4.0', '5.5_3.0']:

    orig = np.loadtxt('d03_%s_A.orig' % model)
    ax.plot(orig[:, 0], orig[:, 1], color='black', lw=1)

    new = np.loadtxt('d03_%s_A.summary' % model)
    ax.plot(new[:, 0], new[:, 2] / new[:, 1], color='red', lw=1)

ax.set_xscale('log')
fig.savefig('comparison_albedo.png')
