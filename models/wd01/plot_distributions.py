import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Loop through R_v values and case A/B
for r_v in [3.1, 4.0, 5.5]:
    for case in ['A', 'B']:

        if r_v == 3.1 and case == 'B':
            continue

        # Open plot and create sub-plots
        fig = plt.figure(figsize=(5, 5))
        ax_s = fig.add_subplot(2, 1, 1)
        ax_g = fig.add_subplot(2, 1, 2)

        # Loop through small carbonaceous grain fractions
        for b_c in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:

            try:

                # Silicate
                s = np.loadtxt('sizes/wd01_silicate_%3.1f_%3.1f_%1s' \
                               % (r_v, b_c, case))
                ax_s.loglog(s[:, 0], s[:, 1] * s[:, 0] ** 4. * 1.e17,\
                            color='black')

                # Graphite
                g = np.loadtxt('sizes/wd01_graphite_%3.1f_%3.1f_%1s' \
                               % (r_v, b_c, case))
                ax_g.loglog(g[:, 0], g[:, 1] * g[:, 0] ** 4. * 1.e17, \
                            color='black')

            except IOError:
                continue

        # Set the limits to resemble those in WD01
        ax_s.set_xlim(1.e-4, 10.)
        ax_s.set_ylim(0.1, 120.)
        ax_g.set_xlim(1.e-4, 10.)
        ax_g.set_ylim(0.1, 120.)

        # Save figure
        fig.savefig('plots/wd01_%3.1f_%1s.png' % (r_v, case))
