import argparse
import json
import math
import numpy as np
import matplotlib.pyplot as plt

from cutoff_estimates import get_max_nbars, get_max_nbars_all_spins, get_cutoff_arr_all_spins, get_cutoff_arr
from timslib.ion_crystals.ion_chain       import IonChain
from timslib.ion_crystals.normal_modes    import nu_ax_from_nu_rad_min
from timslib.ms_pulse_shaping.ms_handler  import MSHandlerPPolyEqF as MSHandler
from fidelity_contribs_from_alpha_and_chi import get_phase_space_results, alpha_spin_string_inf, alpha_avg_inf, delta_chi 
from spin_flip.spin_flip_prob import spin_flip_prob
from fidelity_contribs_from_tdse import get_tdse_results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tdse', action='store_true')
    parser.add_argument('-spin_flip', action='store_true')
    parser.add_argument('-basis', default='x')
    parser.add_argument('-ham_type', default='full')
    args = parser.parse_args()

    n_ions = 2
    nu_rad = 1e6
    nu_rad_min = 0.75e6
    nu_ax = nu_ax_from_nu_rad_min(nu_rad_min, nu_rad, n_ions)

    chain_attrs =  {
        "ion_type" : "Ca40",
        "nu_ax" :  nu_ax,
        "nu_rad" :  nu_rad,
        'n_ions' : n_ions
        }

    n1 = int(math.floor(n_ions/2)) - 1
    n2 = n1 + 1

    shaping_attrs = {
        "anglebypi" : 0.5,
        "mode_type" : "radial",
        "t_gate" : 8.46796804144236e-05,
        "muby2pi" : 1041628.40005989,
        "psibypi" : 0.5,
        "thetabypi" : 0.25,
        "shaping_type" : "3",
        "used_ions" : [n1, n2],
        "n_df" : 2*n_ions + 1,
        "carrier_opt" : False,
    }

    chain = IonChain(**chain_attrs)
    handler = MSHandler.get_from_chain_with_pulse_shaping(chain, **shaping_attrs)
    n_carrier = 1501
    handler.calculate_num_alpha_and_chi_w_carrier(n_carrier)
    phase_space_results = get_phase_space_results(handler)
    print(json.dumps(phase_space_results._asdict(), indent=4))

    s = np.array([1, 1])
    basis = args.basis
    n_tdse = 40001
    ham_type = args.ham_type
    print('max nbars:', get_max_nbars(handler, s))
    print('max nbars (all spins):', get_max_nbars_all_spins(handler))
    #exit()

    if args.spin_flip:
        print('calculating spin flip probability:')
        sf_prob = spin_flip_prob(handler, s) if basis == 'x' else np.nan
        print(f'spin flip probability: {sf_prob}')
    else:
        sf_prob = np.nan

    if basis == 'x':
        inf_full_theory = alpha_spin_string_inf(handler, s) + sf_prob
    elif basis == 'z': 
        inf_full_theory = alpha_avg_inf(handler) + delta_chi(handler)**2

    if args.tdse:
        tol = 1e-8
        n_cutoff_arr = get_cutoff_arr(handler, s, tol)
        print('n_cutoff_arr :', n_cutoff_arr)
        n_cutoff_arr = get_cutoff_arr_all_spins(handler, tol)
        print('n_cutoff_arr (all_spins):', n_cutoff_arr)

        print('solving TDSE:')
        tdse_results = get_tdse_results(handler, s, basis, n_cutoff_arr, n_tdse, ham_type)
        print(json.dumps(tdse_results._asdict(), indent=4))

    print(f'inf_full (theory) : {inf_full_theory}')
