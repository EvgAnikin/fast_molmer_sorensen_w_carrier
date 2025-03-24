from itertools import repeat, product
import json
import math
from multiprocessing import Pool
import numpy as np

from scipy.interpolate import PPoly
from gate_and_pulse_params import get_phase_space_results, PhaseSpaceResultsTuple

from timslib.ion_crystals import IonChain
from timslib.ion_crystals.normal_modes import nu_ax_from_nu_rad_min
from timslib.ms_pulse_shaping.pulse_shaping import ExceedMaxOmegaError
from timslib.ms_pulse_shaping.ms_handler import MSHandlerPPolyEqF as MSHandler


def get_results(chain, shaping_attrs, t_gate, muby2pi, n_vals):
    try:
        handler = MSHandler.get_from_chain_with_pulse_shaping(chain, t_gate=t_gate, muby2pi=muby2pi, **shaping_attrs)
        handler.calculate_num_alpha_and_chi_w_carrier(n_vals)
        res = get_phase_space_results(handler)
    except ExceedMaxOmegaError:
        res = np.nan*np.zeros(len(PhaseSpaceResultsTuple._fields))
    return list(res)


def save_gate_params_mu_t_dep(filename, chain_attrs, shaping_attrs, t_gate_range, muby2pi_range, n_vals, n_workers=10):
    chain = IonChain(**chain_attrs)

    with Pool(n_workers) as p:
        t_prod_range, muby2pi_prod_range = list(zip(*list(product(t_gate_range, muby2pi_range))))
        results = p.starmap(get_results, zip(repeat(chain), repeat(shaping_attrs), 
            t_prod_range, muby2pi_prod_range, repeat(n_vals)))

    n_t = t_gate_range.size
    n_mu = muby2pi_range.size
    n_fields = len(PhaseSpaceResultsTuple._fields)

    result_arr = np.array(results).reshape((n_t, n_mu, n_fields))

    def get_index(key):
        return PhaseSpaceResultsTuple._fields.index(key)

    Omega_max_arr     = result_arr[:, :, get_index('Omega_max')]
    carrier_phase_arr = result_arr[:, :, get_index('carrier_phase')]
    alpha_inf_p_arr   = result_arr[:, :, get_index('alpha_inf_p')]
    alpha_inf_m_arr   = result_arr[:, :, get_index('alpha_inf_m')]
    alpha_inf_avg_arr = result_arr[:, :, get_index('alpha_inf_avg')]
    delta_chi_arr     = result_arr[:, :, get_index('delta_chi')]

    print(filename)
    np.savez(filename,     t_gate_range=t_gate_range,
                           mu_range=2*np.pi*muby2pi_range,
                           Omega_max_arr=Omega_max_arr, 
                           carrier_phase_arr=carrier_phase_arr,
                           alpha_inf_p_arr=alpha_inf_p_arr, 
                           alpha_inf_m_arr=alpha_inf_m_arr, 
                           alpha_inf_avg_arr=alpha_inf_avg_arr, 
                           delta_chi_arr=delta_chi_arr)



if __name__ == '__main__':
    nu_rad = 1e6
    nu_rad_min = 0.75e6

    #For the data used in the manuscript: 
    #dt = 0.1 #(us)
    #delta_mu = 0.001 #(MHz)

    #Rough grid for faster calculation
    dt = 1 #(us)
    delta_mu = 0.01 #(MHz)
    

    muby2pi_range = np.arange(0.5, 1.3, delta_mu)*1e6
    t_gate_range  = np.arange(5, 150, dt)*1e-6

    n_vals = 1001

    shaping_attrs_part = {
            "anglebypi" : 0.5,
            "mode_type" : "radial", 
            "psibypi" : 0.5,
            "thetabypi" : 0.25,
            "shaping_type" : "3",
        }

    n_ions_range = [2] # range(2, 21) to reproduce all data in the manuscript

    for n_ions in n_ions_range:
        for carrier_opt in [False, True]:
            nu_ax = nu_ax_from_nu_rad_min(nu_rad_min, nu_rad, n_ions)

            chain_attrs =  {
                "ion_type" : "Ca40",
                "nu_ax" :  nu_ax,
                "nu_rad" :  nu_rad,
                'n_ions' : n_ions
                }

            n1 = int(math.floor(n_ions/2)) - 1
            n2 = n1 + 1

            shaping_attrs = {**shaping_attrs_part, 
                             'n_df' : 2*n_ions + 1,
                             'used_ions' : [n1,n2], 
                             'carrier_opt' : carrier_opt}

            if carrier_opt:
                filename = f'../data/gate_params_mu_t_dep_opt_{n_ions}ions.npz'
            else:
                filename = f'../data/gate_params_mu_t_dep_noopt_{n_ions}ions.npz'

            n_workers = 4
            save_gate_params_mu_t_dep(filename, chain_attrs, shaping_attrs, t_gate_range, muby2pi_range, n_vals, n_workers)
