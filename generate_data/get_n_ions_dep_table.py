import math
import numpy as np
import pandas as pd
import sqlite3
import sys
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from timslib.ms_pulse_shaping.pulse_shaping import MAX_OMEGABYMU


n_ions_range = np.arange(2,21)

min_points = []
min_points_w_inf = []
inf = []

for n_ions in n_ions_range:
    mu_t_dep_file_noopt = f'../mu_t_dep_data/gate_params_mu_t_dep_noopt_{n_ions}ions.npz'
    mu_t_dep_file_opt = f'../mu_t_dep_data/gate_params_mu_t_dep_opt_{n_ions}ions.npz'

    datafile     = np.load(mu_t_dep_file_noopt)
    datafile_opt = np.load(mu_t_dep_file_opt)
    mu_range = datafile['mu_range']
    t_gate_range = datafile['t_gate_range']
    alpha_avg_inf_arr = datafile_opt['alpha_inf_avg_arr']
    delta_chi_arr = datafile_opt['delta_chi_arr']
    alpha_avg_inf_arr_noopt = datafile['alpha_inf_avg_arr']
    delta_chi_arr_noopt = datafile['delta_chi_arr']
    Omega_max_arr = datafile['Omega_max_arr']
    inf_arr = delta_chi_arr**2 + alpha_avg_inf_arr
    allowed_points = []
    allowed_points_inf = []
    
    for i in range(len(t_gate_range)):
        for j in range(len(mu_range)):
            if Omega_max_arr[i][j]/mu_range[j] <= MAX_OMEGABYMU and mu_range[j] > 2*np.pi*1e6:
                allowed_points.append([t_gate_range[i], mu_range[j]])
                if alpha_avg_inf_arr[i][j] + delta_chi_arr[i][j]**2 < 1e-5:
                    allowed_points_inf.append([t_gate_range[i], mu_range[j]])
                    
    allowed_points.sort(key = lambda x: x[0])
    allowed_points_inf.sort(key = lambda x: x[0])
    min_points.append(allowed_points[0])
    min_points_w_inf.append(allowed_points_inf[0])
    j = np.where(mu_range==allowed_points_inf[0][1])[0][0]
    i = np.where(t_gate_range==allowed_points_inf[0][0])[0][0]
    inf.append(delta_chi_arr_noopt[i][j]**2 + alpha_avg_inf_arr_noopt[i][j])

min_t_gate_arr = [pt[0] for pt in min_points]
min_mu_arr     = [pt[1] for pt in min_points]

min_t_gate_w_inf = [pt[0] for pt in min_points_w_inf]
min_mu_w_inf = [pt[1] for pt in min_points_w_inf]

noopt_inf = inf

df = pd.DataFrame({'n_ions' : n_ions_range, 'min_t_gate' : min_t_gate_arr, 'min_mu' : min_mu_arr, 'min_t_gate_1e-5': min_t_gate_w_inf, 'min_mu_1e-5' : min_mu_w_inf, 'noopt_inf' : noopt_inf})

print(df)
    
#ax = plt.figure(figsize=(5,3)).gca()
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#plt.ylabel(r'gate time ($\mu$s)',fontsize=9)
#plt.xlabel(r'$n_{ions}$',fontsize=9)
#plt.plot(n_ions, [1e6*pt[0] for pt in min_points],'o' ,label = r'$t_{min}$')
#plt.plot(N_ions, [1e6*pt[0] for pt in min_points_w_inf],'>' ,label = r'$t_{min}^{*}$')
#plt.legend(fontsize=9)
#plt.savefig('Figure_1.png', bbox_inches='tight')
#
#
#ax1 = plt.figure(figsize=(5,3)).gca()
#ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
#plt.xlabel(r'$N_{ions}$',fontsize=9)
#plt.ylabel(r'$1-F_{0}$',fontsize=9)
#plt.ylim(bottom=1e-5)
#plt.plot(N_ions, inf,'o' ,label = 'noopt inf')
#plt.yscale('log')
#plt.legend(fontsize=9)
#plt.savefig('Figure_2.png', bbox_inches='tight')
#plt.show()
#plt.show()

#print([pt[1] for pt in min_points_w_inf])
