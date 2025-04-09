import argparse
import collections
import datetime
import itertools
import math
import matplotlib
import matplotlib.gridspec as gridspec
#import matplotlib.colormaps
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sqlite3
import timeit
import sys

from scipy.interpolate import CubicSpline, PPoly


def plot_ms_inf_z(ax):
    db_name = '../data/ms_data.db'

    with sqlite3.connect(db_name) as db:
        query= """SELECT t_gate, muby2pi, inf_theory, inf_num FROM 
                                (MU_T_POINTS INNER JOIN RESULTS USING(point_id)) 
                                WHERE n_ions == {}
                                AND basis == 'z'
                                AND init_spin_2 == 1 
                                AND ham_type == {}
                                AND carrier_opt == {}"""


        df_noopt_ld   = pd.read_sql_query(query.format(5, "'LD'",   "FALSE"), db)
        df_opt_ld     = pd.read_sql_query(query.format(5, "'LD'", "TRUE"),  db)
        df_noopt_full = pd.read_sql_query(query.format(5, "'full'",   "FALSE"), db)
        df_opt_full   = pd.read_sql_query(query.format(5, "'full'", "TRUE"),  db)


    ax.scatter(df_noopt_ld['t_gate']*1e6,   df_noopt_ld['inf_theory'], label='$1-F_0$ (lin)', marker='>')
    ax.scatter(df_opt_ld['t_gate']*1e6,     df_opt_ld['inf_theory'],   label='$1-F_0$ (tr)', marker='<')
    ax.scatter(df_noopt_ld['t_gate']*1e6,   df_noopt_ld['inf_num'], label='lin, LD', marker='x')
    ax.scatter(df_opt_ld['t_gate']*1e6,     df_opt_ld['inf_num'], label='tr, LD', marker='+')
    ax.scatter(df_noopt_full['t_gate']*1e6, df_noopt_full['inf_num'], label='lin, full', marker='*')
    ax.scatter(df_opt_full['t_gate']*1e6,   df_opt_full['inf_num'], label='tr, full', marker='h')


    ax.legend(loc='right')
    ax.set_yscale('log')
    ax.set_xlim(38, 150)


if __name__ == '__main__':
    fig, ax = plt.subplots(1,1, figsize=(5,4))

    ax.set_xlabel('$t_{gate}$ (us)')
    ax.set_ylabel('$1-F_0$')

    plot_ms_inf_z(ax)

    plt.savefig('../figures/z_states_inf.pdf')
    plt.tight_layout()
    plt.show()
