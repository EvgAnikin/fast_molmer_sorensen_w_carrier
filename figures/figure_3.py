import argparse
import collections
import datetime
import itertools
import math
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sqlite3
import timeit
import sys
from matplotlib.colors import LogNorm
from scipy.interpolate import CubicSpline, PPoly

from timslib.ion_crystals import IonChain
from timslib.ion_crystals.normal_modes import nu_ax_from_nu_rad_min

MAX_OMEGABYMU = 0.5818652242815963

def indicate_normal_modes(ax, nus, x0, dx, **kwargs):
    for nu in nus:
        ax.plot([x0, x0+dx], [nu, nu], **kwargs)

def plot_value_mu_t_dep_on_ax(ax, n_ions, value_arr_name, carrier_opt, norm, cmap):
    name_noopt = f'../data/gate_params_mu_t_dep_noopt_{n_ions}ions.npz'
    name_opt   = f'../data/gate_params_mu_t_dep_opt_{n_ions}ions.npz'

    if carrier_opt:
        datafile = np.load(name_opt)
        datafile_noopt = np.load(name_noopt)
    else:
        datafile = np.load(name_noopt)
        datafile_noopt = datafile

    mu_range = datafile['mu_range']
    t_gate_range = datafile['t_gate_range']
    Omega_max_arr = datafile_noopt['Omega_max_arr']
    alpha_inf_avg_arr = datafile['alpha_inf_avg_arr']
    delta_chi_arr = datafile['delta_chi_arr']

    if value_arr_name == 'alpha_inf_avg_arr':
        value_arr = alpha_inf_avg_arr
    elif value_arr_name == 'delta_chi_arr':
        value_arr = delta_chi_arr
    elif value_arr_name == 'inf_z':
        value_arr = delta_chi_arr**2 + alpha_inf_avg_arr

    if value_arr_name == 'Omega_max_arr':
        ax.contour(t_gate_range*1e6, mu_range/(2*math.pi)*1e-6, (Omega_max_arr/mu_range[np.newaxis, :]).T, levels=[MAX_OMEGABYMU], colors=['k'])
        ax.contourf(t_gate_range*1e6, mu_range/(2*math.pi)*1e-6, (Omega_max_arr/mu_range[np.newaxis, :]).T, levels=[0, MAX_OMEGABYMU], colors=['k'], alpha=0.4)
    else:
        CS = ax.pcolormesh(t_gate_range*1e6, mu_range/(2*math.pi)*1e-6, np.abs(value_arr.T), norm=norm, cmap=cmap)

    nu_rad = 1e6
    nu_rad_min = 0.75e6
    nu_ax = nu_ax_from_nu_rad_min(nu_rad_min, nu_rad, n_ions)
    chain = IonChain(n_ions=n_ions, nu_ax=nu_ax, nu_rad=nu_rad, ion_type='Ca40')
    indicate_normal_modes(ax, chain.omegas_rad/(2*np.pi)/1e6, t_gate_range[0]*1e6, 5, color='k', linewidth=0.8)


def plot_mu_t_points_on_ax(ax, n_ions, **scatter_kwargs):
    db_name = '../data/ms_data.db'
    with sqlite3.connect(db_name) as db:
        query = """SELECT t_gate, muby2pi, spin_flip_theory FROM 
                                (MU_T_POINTS INNER JOIN RESULTS USING(point_id)) 
                                WHERE n_ions == {} AND muby2pi > 1e6"""
        df = pd.read_sql_query(query.format(n_ions), db)
    ax.scatter(df['t_gate']*1e6, df['muby2pi']*1e-6, **scatter_kwargs)

def plot_carrier_error_5ions(ax):
    db_name = '../data/ms_data.db'
    with sqlite3.connect(db_name) as db:
        query= """SELECT t_gate, muby2pi, spin_flip_theory, spin_flip_num FROM 
                                (MU_T_POINTS INNER JOIN RESULTS USING(point_id)) 
                                WHERE n_ions == {}
                                AND basis == 'x'
                                AND init_spin_2 == 1 
                                AND ham_type == {}
                                AND carrier_opt == {}"""
        df_noopt_ld   = pd.read_sql_query(query.format(5, "'LD'",   "FALSE"), db)
        df_opt_ld     = pd.read_sql_query(query.format(5, "'LD'", "TRUE"),  db)
        df_noopt_full = pd.read_sql_query(query.format(5, "'full'",   "FALSE"), db)
        df_opt_full   = pd.read_sql_query(query.format(5, "'full'", "TRUE"),  db)
    ax.scatter(df_noopt_ld['t_gate']*1e6, df_noopt_ld['spin_flip_theory'], label=r'$\Omega_\mathrm{lin}$', marker='>')
    ax.scatter(df_opt_ld['t_gate']*1e6, df_opt_ld['spin_flip_theory'], 
            label=r'$\Omega_\mathrm{tr}$', marker='<')
    ax.scatter(df_noopt_ld['t_gate']*1e6, df_noopt_ld['spin_flip_num'], 
            label=r'$\Omega_\mathrm{lin}$, LD', marker='x')
    ax.scatter(df_opt_ld['t_gate']*1e6,   df_opt_ld['spin_flip_num'],     label=r'$\Omega_\mathrm{tr}$, LD', marker='+')
    ax.scatter(df_noopt_ld['t_gate']*1e6, df_noopt_full['spin_flip_num'], label=r'$\Omega_\mathrm{lin}$, full',         marker='*')
    ax.scatter(df_opt_ld['t_gate']*1e6,   df_opt_full['spin_flip_num'],   label=r'$\Omega_\mathrm{tr}$, full',           marker='h')
    ax.legend(loc='lower right')
    ax.set_yscale('log')

def plot_carrier_error_20ions(ax):
    db_name = '../data/ms_data.db'
    with sqlite3.connect(db_name) as db:
        query = """SELECT t_gate, muby2pi, spin_flip_theory FROM 
                                (MU_T_POINTS INNER JOIN RESULTS USING(point_id)) 
                                WHERE n_ions == {}
                                AND muby2pi > 1e6
                                AND basis == 'x'
                                AND init_spin_2 == 1 
                                AND ham_type == 'LD'
                                AND carrier_opt == {}"""
        df_noopt = pd.read_sql_query(query.format(20, "FALSE"), db)
        df_opt = pd.read_sql_query(query.format(20, "TRUE"), db)
    ax.scatter(df_noopt['t_gate']*1e6, df_noopt['spin_flip_theory'], 
            label=r'$\Omega_\mathrm{lin}$', marker='>')
    ax.scatter(df_opt['t_gate']*1e6, df_opt['spin_flip_theory'], 
            label=r'$\Omega_\mathrm{tr}$', marker='<')
    ax.legend()
    ax.set_yscale('log')

if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 11})
    fig = plt.figure(figsize=(10, 10)) 

    gs = gridspec.GridSpec(3, 2)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[1, 0])
    ax4 = plt.subplot(gs[1, 1])
    ax5 = plt.subplot(gs[2, 0])
    ax6 = plt.subplot(gs[2, 1])

    norm = LogNorm(vmin=1e-6, vmax=1e-1)
    cmap = plt.get_cmap('gnuplot')
    
    plot_value_mu_t_dep_on_ax(ax1, 5, 'inf_z',  carrier_opt=True,  norm=norm, cmap=cmap)
    plot_value_mu_t_dep_on_ax(ax3, 5, 'inf_z',  carrier_opt=False, norm=norm, cmap=cmap)
    plot_value_mu_t_dep_on_ax(ax2, 20, 'inf_z', carrier_opt=True,  norm=norm, cmap=cmap)
    plot_value_mu_t_dep_on_ax(ax4, 20, 'inf_z', carrier_opt=False, norm=norm, cmap=cmap)

    plot_mu_t_points_on_ax(ax1, 5,  color='red', zorder=3, s=10, marker='D')
    plot_mu_t_points_on_ax(ax2, 20, color='red', zorder=3, s=10, marker='D')

    plot_carrier_error_5ions(ax5)
    plot_carrier_error_20ions(ax6)

    ax5.set_xlabel('$t_{gate}$ ($\mu$s)')
    ax6.set_xlabel('$t_{gate}$ ($\mu$s)')
    ax1.set_ylabel('$\mu/(2\pi)$ (MHz)')
    ax3.set_ylabel('$\mu/(2\pi)$ (MHz)')
    ax5.set_ylabel('$\mu/(2\pi)$ (MHz)')
    ax5.set_xlim(5,150)
    ax6.set_xlim(5,150)
    ax5.set_ylim(1e-9,1e-4)
    ax6.set_ylim(1e-9,1e-4)

    labels = [f'({c})' for c in 'abcdef']
    for ax, label in zip([ax1, ax2, ax3, ax4, ax5, ax6], labels):
        ax.text(0.05, 0.9, label, transform=ax.transAxes)
    
    cbar_ax = fig.add_axes([0.15, 0.96, 0.7, 0.017])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar_ax.set_xlabel("$1-F_0$")
    cbar_ax.xaxis.set_label_position('top')
    
    plt.tight_layout(rect=[0, 0, 1, 0.94])  
    plt.show()
