from collections import namedtuple
import numpy as np
from ions_tdse_qutip import chain_spin_string, get_tdse_solution
from qutip import *


def state_from_arr(state_arr, n_ions, n_cutoff_arr):
    n_modes = n_cutoff_arr.size
    state = Qobj(state_arr.T[0], dims=([n_ions*[2] + list(n_cutoff_arr), (n_ions + n_modes)*[1]]))
    return state


def phonon_excitation_no_spin_flip(psi_f, s, n_ions, n_cutoff_arr, basis):
    if len(s) > 2:
        raise ValueError
    
    psi_0_spin = chain_spin_string(s, [], spin_basis=basis)
    psi_0      = chain_spin_string(s, n_cutoff_arr, spin_basis=basis)   
    rho = psi_f.ptrace([0,1])
    
    P0  = np.real((psi_0_spin.dag()*rho*psi_0_spin))
    P00 = np.abs((psi_0.dag()*psi_f))**2
        
    return P0 - P00


def spin_flip_excitation(psi, s, basis):
    rho = psi.ptrace([0,1])


    if len(s) > 2:
        raise ValueError
    
    s1 = np.array([ s[0], -s[1]])
    s2 = np.array([-s[0],  s[1]])
 
    psi_1 = chain_spin_string(s1, [], spin_basis=basis)
    psi_2 = chain_spin_string(s2, [], spin_basis=basis)

    P1 = (psi_1.dag()*rho*psi_1).real
    P2 = (psi_2.dag()*rho*psi_2).real
    
    print(f'P1 : {P1}')
    print(f'P2 : {P2}')
    
    return P1 + P2


def MS_full_fidelity(real_state, psi_0_q, n_ions, n_cutoff_arr):
    ideal_MS_p = ( 1j*np.pi/4*tensor(sigmax(), sigmax())).expm()
    ideal_MS_m = (-1j*np.pi/4*tensor(sigmax(), sigmax())).expm()
    ideal_state_p = tensor([ideal_MS_p * psi_0_q] + [fock(n_cutoff, 0) for n_cutoff in n_cutoff_arr])
    ideal_state_m = tensor([ideal_MS_m * psi_0_q] + [fock(n_cutoff, 0) for n_cutoff in n_cutoff_arr])
    return max(np.abs((ideal_state_p.dag()*real_state))**2, np.abs((ideal_state_m.dag()*real_state))**2)


def acquired_phase(final_state, s, basis, n_ions, n_cutoff_arr):
    psi_0_full = chain_spin_string(s, n_cutoff_arr, spin_basis=basis)
    return np.angle((psi_0_full.dag()*final_state))


TDSEResultTuple = namedtuple('TDSEResultTuple', 
        ['fin_phase_num', 'ph_ex_num', 'spin_flip_num', 'inf_num'])

def get_tdse_results(handler, s, basis, n_cutoff_arr, n_t, ham_type='full'):
    t_range, states = get_tdse_solution(handler, s, basis, n_cutoff_arr, n_t, ham_type)
    fin_state = states[-1]

    n_ions = handler.n_used_ions
    psi0_q = chain_spin_string(s, [], spin_basis=basis)

    inf_full = 1 - MS_full_fidelity(fin_state, psi0_q, n_ions, n_cutoff_arr)

    if basis == 'x':
        ph_ex_num     = phonon_excitation_no_spin_flip(fin_state, s, n_ions, n_cutoff_arr, basis)
        sf_ex_num     = spin_flip_excitation(fin_state, s, basis)
        fin_phase_num = acquired_phase(fin_state, s, basis, n_ions, n_cutoff_arr)
    else:
        ph_ex_num = np.nan
        sf_ex_num = np.nan
        fin_phase_num = np.nan

    return TDSEResultTuple(fin_phase_num, ph_ex_num, sf_ex_num, inf_full)
