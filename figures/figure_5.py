import collections
import datetime
import itertools
import math
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import timeit
import sys

from matplotlib.ticker import MaxNLocator



if __name__ == '__main__':
    df = pd.read_csv('../data/n_ions_dep.txt', sep='\s+')

    matplotlib.rcParams.update({'font.size': 11})
#    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(5,6))
    fig, ax1 = plt.subplots(1,1, figsize=(5,3.5))

#    ax1.text(0.04, 0.9, '(a)', transform=ax1.transAxes, fontsize=12)
    ax1.scatter(df['n_ions'], df['min_t_gate'], marker='x', label='$t_{min}$')
    ax1.scatter(df['n_ions'], df['min_t_gate_1e-5'], marker='>', label='$t_{min}^*$')
    ax1.set_xlabel('$n_{ions}$')
    ax1.set_ylabel('gate time ($\mu$s)')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.legend()

#    ax2.text(0.04, 0.9, '(b)', transform=ax2.transAxes, fontsize=12)
#    ax2.scatter(df['n_ions'], df['noopt_inf'])
#    ax2.set_ylim(1e-3, 1e-1)
#    ax2.set_yscale('log')
#    ax2.set_xlabel('$n_{ions}$')
#    ax2.set_ylabel('$1-F_0$')
#    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig('../figures/n_ions_dep.pdf')
    plt.show()
