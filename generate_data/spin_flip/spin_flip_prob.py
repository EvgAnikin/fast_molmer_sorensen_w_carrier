import numpy as np
from . import spin_flip


def spin_flip_prob(handler, s):
    n_vals =  handler.alpha_num_w_carrier.shape[0]
    t_range = np.linspace(handler.t_i, handler.t_f, n_vals)
    Omega_range = handler.Omega(t_range)
    carrier_phases = np.array([handler.carrier_phase(t, handler.t_i) for t in t_range])
    return spin_flip.spin_flip_prob(s, t_range, Omega_range, handler.eta, handler.omegas, handler.mu, handler.psi, 
            carrier_phases, handler.mode_alpha_num_w_carrier, handler.mode_chi_num_w_carrier)
