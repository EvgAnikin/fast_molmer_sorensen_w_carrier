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
