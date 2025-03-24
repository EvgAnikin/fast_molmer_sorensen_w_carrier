import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy

from qutip import *


def local_operator(op, m, dimensions):
    if op.shape[0] != op.shape[1] or op.shape[0] != dimensions[m]:
        raise ValueError('Incompatible dimensions')

    operators_list = [qeye(dim) for dim in dimensions]
    operators_list[m] = op

    return tensor(operators_list)


def chain_gs(n_cutoff_arr):
    return tensor(fock(n_cutoff, 0) for n_cutoff in n_cutoff_arr)

def spin_eigenstate(s, spin_basis):
    if spin_basis == 'z':
        if s == 1:
            return fock(2,0)
        elif s == -1:
            return fock(2,1)
    elif spin_basis == 'x':
        if s == 1:
            return 1/math.sqrt(2)*(fock(2, 0) + fock(2, 1))
        elif s == -1:
            return 1/math.sqrt(2)*(fock(2, 0) - fock(2, 1))
    elif spin_basis == 'y':
        if s == 1:
            return 1/math.sqrt(2)*(fock(2, 0) + 1j*fock(2, 1))
        elif s == -1:
            return 1/math.sqrt(2)*(1j*fock(2, 0) + fock(2, 1))


def chain_spin_string(s, n_cutoff_arr, spin_basis='z'):
    ion_states_list = [spin_eigenstate(s_i, spin_basis) for s_i in s]
    phonon_states_list = [fock(n_cutoff, 0) for n_cutoff in n_cutoff_arr]

    return tensor(ion_states_list + phonon_states_list)


def chain_ion_op(operator, i, n_ions, n_cutoff_arr):
    dimensions = np.concatenate((2*np.ones(n_ions, dtype=int), n_cutoff_arr)).tolist()
    return local_operator(operator, i, dimensions)


def chain_sigmax(i, n_ions, n_cutoff_arr):
    return chain_ion_op(sigmax(), i, n_ions, n_cutoff_arr)


def chain_sigmay(i, n_ions, n_cutoff_arr):
    return chain_ion_op(sigmay(), i, n_ions, n_cutoff_arr)


def chain_sigmaz(i, n_ions, n_cutoff_arr):
    return chain_ion_op(sigmaz(), i, n_ions, n_cutoff_arr)


def chain_sigmap(i, n_ions, n_cutoff_arr):
    return chain_ion_op(sigmap(), i, n_ions, n_cutoff_arr)


def chain_sigmam(i, n_ions, n_cutoff_arr):
    return chain_ion_op(sigmam(), i, n_ions, n_cutoff_arr)


def chain_phonon_op(operator, m, n_ions, n_cutoff_arr):
    dimensions = np.concatenate((2*np.ones(n_ions, dtype=int), n_cutoff_arr)).tolist()
    return local_operator(operator, n_ions + m, dimensions)


def chain_create(m, n_ions, n_cutoff_arr):
    return chain_phonon_op(create(n_cutoff_arr[m]), m, n_ions, n_cutoff_arr)


def chain_create(m, n_ions, n_cutoff_arr):
    return chain_phonon_op(destroy(n_cutoff_arr[m]), m, n_ions, n_cutoff_arr)


def chain_q(m, n_ions, n_cutoff_arr):
    q = create(n_cutoff_arr[m]) + destroy(n_cutoff_arr[m])
    return chain_phonon_op(q, m, n_ions, n_cutoff_arr)


def chain_num(m, n_ions, n_cutoff_arr):
    return chain_phonon_op(num(n_cutoff_arr[m]), m, n_ions, n_cutoff_arr)


def exp_unitary(op):
    energies, vectors = scipy.linalg.eigh(1j*op)
    eigvals = np.exp(-1j*energies)
    return np.matmul(vectors, np.matmul(np.diag(eigvals), np.conj(vectors.T)))


def ion_ham(handler, n_cutoff_arr, n_points_t, ham_type=False):
    if handler.n_modes != len(n_cutoff_arr): 
        raise ValueError('Invalid n_cutoff_arr size')

    free_phonon_ham = sum([handler.omegas[m]*chain_num(m, handler.n_used_ions, n_cutoff_arr) for m in range(handler.n_modes)])

    t_range = np.linspace(handler.t_i, handler.t_f, n_points_t)
    Omegas = handler.Omega(t_range)
     

    n_ions = handler.n_used_ions
    n_modes = handler.n_modes

    q_ops = [chain_q(m, n_ions, n_cutoff_arr) for m in range(len(n_cutoff_arr))]

    ion_phonon_ham = 0
    for i in range(n_ions):
        eta = handler.eta
        eta_q = sum([eta[i,m]*q_ops[m] for m in range(n_modes)])

        if ham_type == 'MS':
            exp_ietaq = 1j*eta_q
        elif ham_type == 'LD':
            exp_ietaq = 1 + 1j*eta_q 
        elif ham_type == 'full':
            exp_ietaq = Qobj(exp_unitary(1j*eta_q.data.to_array()), dims=eta_q.dims)
        else:
            raise ValueError('Unknown ham_type: should be either MS, LD or full')

        cos_etaq = Qobj( 0.5 *(exp_ietaq + exp_ietaq.dag()))
        sin_etaq = Qobj(-0.5j*(exp_ietaq - exp_ietaq.dag()))

        sigmax_i = chain_sigmax(i, n_ions, n_cutoff_arr)
        sigmay_i = chain_sigmay(i, n_ions, n_cutoff_arr)

        if ham_type == 'pure_MS':
            ion_phonon_ham += sigmax_i*eta_q
        else:
            ion_phonon_ham += sigmay_i*cos_etaq + sigmax_i*sin_etaq 

    ham = [free_phonon_ham, [ion_phonon_ham,  Omegas*np.cos(handler.mu*t_range + handler.psi)]]
    return QobjEvo(ham, tlist=t_range)


def get_tdse_solution(handler, s, basis, n_cutoff_arr, n_t, ham_type=False):
    ham = ion_ham(handler, n_cutoff_arr=n_cutoff_arr, n_points_t=n_t, ham_type=ham_type)
    t_range = np.linspace(handler.t_i, handler.t_f, n_t)
    n_ions = handler.n_used_ions
    psi_0 = chain_spin_string(s, n_cutoff_arr, basis)
    result = sesolve(ham, psi_0, tlist=(handler.t_i, handler.t_f), options={'nsteps': n_t, 'rtol' : 1e-12, 'atol' : 1e-12})
    states = result.states
    
    return t_range, states
