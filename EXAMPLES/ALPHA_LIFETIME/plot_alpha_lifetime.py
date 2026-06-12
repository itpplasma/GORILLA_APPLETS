#!/usr/bin/env python3
"""
Plot confined alpha-particle fraction vs. tracing time.

Reads alpha_lifetime_gorilla.dat (one confinement time per line, in seconds)
produced by a prior run of gorilla_applets_main.x with i_option = 5.

Python translation of plot_alpha_lifetime.m.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcdefaults()


t_confined = np.atleast_1d(np.genfromtxt('alpha_lifetime_gorilla.dat'))
n_sampling = t_confined.size

t_sample = np.logspace(-4, np.log10(t_confined.max()), 100)
confined_ratio = np.array([(t_confined >= ts).sum() / n_sampling for ts in t_sample])
confined_sigma = np.sqrt(confined_ratio * (1.0 - confined_ratio) / n_sampling)

fig, ax = plt.subplots(figsize=(10, 8))
ax.fill_between(t_sample,
                confined_ratio - 1.96 * confined_sigma,
                confined_ratio + 1.96 * confined_sigma,
                color='k', alpha=0.2, linewidth=0)
ax.plot(t_sample, confined_ratio, color='r', linewidth=2,
        label=r'GORILLA Poly4, $N_s = N_\vartheta = N_\varphi = 101$')

ax.set_xscale('log')
ax.set_xlabel(r'$t_\mathrm{trace}$ [s]')
ax.set_ylabel(r'$f_c$')
ax.set_ylim(0.9, 1.0)
ax.grid(True, which='both')
ax.set_title(r'$t_\mathrm{trace} = 0.01$ s, $N_\mathrm{particles} = 10^4$')
ax.legend(loc='lower left')

plt.tight_layout()
plt.savefig('alpha_lifetime.png', dpi=150)
plt.show()
