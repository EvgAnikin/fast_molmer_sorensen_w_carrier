from collections import namedtuple
import itertools
import math
import numpy as np
from scipy.interpolate import PPoly

def Omega_max_abs(t_bp, Omega_mat):
    def spline_extrema(spline):
        return spline.derivative().roots()
    
    def spline_max_abs(spline):
        return np.max(np.abs(spline(spline_extrema(spline))))

    Omega = PPoly.construct_fast(Omega_mat, t_bp)
    return spline_max_abs(Omega)


def alpha_spin_string_inf(handler, s):
    alpha_s = np.dot(handler.alpha_fin_num_w_carrier.T, s)
    return np.sum(np.abs(alpha_s)**2)


def alpha_avg_inf(handler):
    return np.sum(np.abs(handler.alpha_fin_num_w_carrier)**2)


def delta_chi(handler):
    def closest_pm_delta(a, b):
        return min(abs(a-b), abs(a+b))
    
    return closest_pm_delta(handler.chi_fin_num_w_carrier[0,1], np.pi/4)


def get_max_nbars(handler, s):
    alpha = np.tensordot(handler.alpha_num_w_carrier, s, axes=[1,0])
    nbars = np.max(np.abs(alpha)**2, axis=0)
    return nbars


def get_max_nbars_all_spins(handler):
    nbars_arr = []

    for s_tuple in itertools.product([1,-1], repeat=handler.n_used_ions):
        s = np.array(s_tuple)
        nbars_arr.append(get_max_nbars(handler, s))

    nbars_arr = np.array(nbars_arr)
    return np.max(nbars_arr, axis=0)


@np.vectorize
def cutoff_from_nbar(nbar, tol, min_cutoff, max_cutoff):
    def p(alpha, n):
        return math.exp(-nbar)*nbar**n/math.factorial(n)
    n = round(nbar)
    while p(nbar, n) > tol:
        n += 1

    if n < min_cutoff:
        return min_cutoff
    elif n > max_cutoff:
        return max_cutoff
    else:
        return n


def get_cutoff_arr(handler, s, tol, min_cutoff=1, max_cutoff=100):
    nbars = get_max_nbars(handler, s)
    return cutoff_from_nbar(nbars, tol, min_cutoff, max_cutoff)


def get_cutoff_arr_all_spins(handler, tol, min_cutoff=1, max_cutoff=100):
    nbars = get_max_nbars_all_spins(handler)
    return cutoff_from_nbar(nbars, tol, min_cutoff, max_cutoff)


PhaseSpaceResultsTuple = namedtuple('PhaseSpaceResultsTuple', 
        ['Omega_max', 'carrier_phase', 'alpha_inf_p', 'alpha_inf_m', 'alpha_inf_avg', 'delta_chi'])

def get_phase_space_results(handler):

    Omega_max = Omega_max_abs(handler.t_bp, handler.Omega_mat)
    carrier_phase = handler.carrier_phases_at_bp[-1]
    alpha_inf_p = alpha_spin_string_inf(handler, np.array([1,1]))
    alpha_inf_m = alpha_spin_string_inf(handler, np.array([1,-1]))
    alpha_inf_avg = alpha_avg_inf(handler)
    dchi = delta_chi(handler)

    return PhaseSpaceResultsTuple(Omega_max, carrier_phase, alpha_inf_p, alpha_inf_m, alpha_inf_avg, dchi)
