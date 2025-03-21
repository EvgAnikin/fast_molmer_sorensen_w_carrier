import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


from .seg_ppoly_transform import *
from .ms_evo_funcs.ppoly_seg_matrices import *

from scipy.special import j0, jv
from scipy.optimize import fsolve


def matrix_form_min_power_segmentation(theta, omegas, alpha_matrix, xi_matrix, norm_matrix, inexact_closing_allowed=True):
    alpha_re_im = np.block([[alpha_matrix.real], [alpha_matrix.imag]])  
    null_vectors = sp.linalg.null_space(alpha_re_im)

    if null_vectors.size == 0: 
        if not inexact_closing_allowed:
            raise ValueError('null subspace is empty: not enough segments/degrees of freedom. alpha_re_im shape: {}'.format(alpha_re_im.shape))
        alpha_err_matrix = np.real(np.conj(alpha_matrix.T) @ np.diag(1/omegas) @ alpha_matrix)
        eigvals, vectors = sp.linalg.eigh(alpha_err_matrix, norm_matrix)
        min_vector = vectors[:, 0]
        multiplier = min_vector.T @ xi_matrix @ min_vector
        
        if multiplier >= 0:
            result = min_vector*math.sqrt(theta)/math.sqrt(multiplier)
        else:
            result = min_vector*math.sqrt(np.pi/2 - theta)/math.sqrt(-multiplier)
        return result
        

    Y_on_ns = null_vectors.T @ xi_matrix  @ null_vectors
    N_on_ns = null_vectors.T @ norm_matrix @ null_vectors

    eigvals, vectors = sp.linalg.eigh(Y_on_ns, N_on_ns)

    min_eigval = eigvals[0]
    min_vector = vectors[:, 0]

    max_eigval = eigvals[-1]
    max_vector = vectors[:, -1]

    which_eigval = 'max'
    if min_eigval >= 0:
        which_eigval = 'max'
    elif max_eigval <= 0:
        which_eigval = 'min'
    elif min_eigval < 0 and max_eigval > 0:
#        if max_eigval/abs(min_eigval) > theta/(math.pi/2-theta):
        if max_eigval*(math.pi/2 - theta) > theta*abs(min_eigval):
            which_eigval = 'max'
        else:
            which_eigval = 'min'

    if which_eigval == 'max':
        eigval = max_eigval
        vector = max_vector
        theta_res = theta
    elif which_eigval == 'min':
        eigval = min_eigval
        vector = min_vector
        theta_res = math.pi/2 - theta

    result = math.sqrt(theta_res)/math.sqrt(np.abs(eigval))*null_vectors @ vector
    return result


def get_ppoly_pulse_from_t_bp(theta, t_bp, eta1, eta2, mu, psi, omegas, kind='pointwise spline'):

    if kind in ['pointwise spline', '3']:
        T = pointwise_spline_matrix(t_bp)
    elif kind in ['smoothed stepwise', '0s']:
        T = smoothed_piecewise_matrix(t_bp)
    elif kind in ['linear', '1']:
        T = linear_interp_matrix(t_bp)
    elif kind in ['constant piecewise', '0']:
        T = constant_piecewise_matrix(t_bp)
    else:
        raise ValueError("""Unknown segmentation type {}: must be either \'pointwise spline\' or \'3\', 
                            \'smoothed stepwise\' or \'0s\', \'linear\' or \'1\', \'constant piecewise\' or \'0\'""")

    # The matrix T_{psk}, where p - polynomial index, s - segment number, k - basis function number

    pmax = T.shape[0] - 1
    n_seg = t_bp.size - 1
    N = ppoly_norm_matrix_reversed(t_bp, pmax)

    alpha_sm = alpha_seg_matrix(t_bp, mu, psi, omegas, pmax) # The matrix A_{mps} (m - mode number, p - polynomial index, s - segment number)
    x_sm     = x_seg_matrix(t_bp, eta1, eta2, mu, psi, omegas, pmax)

    alpha_new_basis = np.tensordot(alpha_sm, T, 2)
    x_new_basis = np.tensordot(np.tensordot(T, x_sm, axes=[[0, 1], [0, 1]]), T)
    N_new_basis = np.tensordot(np.tensordot(T, N, axes=[[0, 1], [0, 1]]), T)

    segmentation_vector = matrix_form_min_power_segmentation(theta, omegas, alpha_new_basis, x_new_basis, N_new_basis)
    Omega_mat = T @ segmentation_vector
    return Omega_mat


def get_ppoly_pulse(eta1, eta2, omegas, muby2pi, t_gate, psibypi, thetabypi, n_df, shaping_type, smoothing_degree=0.1):
    mu = 2*math.pi*muby2pi
    psi = math.pi*psibypi
    theta = math.pi*thetabypi


    if shaping_type in ['pointwise spline', '3']:
        t_bp = np.linspace(0, t_gate, n_df + 2)
    elif shaping_type in ['smoothed stepwise', '0s']:
        p = smoothing_degree
        t_plat =  (1 - p)/(n_df + p)*t_gate
        t_trans = p/(n_df + p)*t_gate
        t_bp = breakpoints_for_smoothed_piecewise(n_df, t_plat, t_trans)
    elif shaping_type in ['linear', '1']:
        t_bp = np.linspace(0, t_gate, n_df + 2)
    elif shaping_type in ['constant piecewise', '0']:
        t_bp = np.linspace(0, t_gate, n_df + 1)
    else:
        raise ValueError("""Unknown segmentation type {}: must be either \'pointwise spline\' or \'3\', 
                            \'smoothed stepwise\' or \'0s\', \'linear\' or \'1\', \'constant piecewise\' or \'0\'""")

    return t_bp, get_ppoly_pulse_from_t_bp(theta, t_bp, eta1, eta2, mu, psi, omegas, shaping_type)


def get_ppoly_pulse_from_chain(chain, anglebypi, mode_type, used_ions, *args, **kwargs):
    angle = math.pi*anglebypi

    if mode_type == 'axial':
        eta = chain.eta_ax(angle)
        omegas = chain.omegas_ax
    elif mode_type == 'radial':
        eta = chain.eta_rad(angle)
        omegas = chain.omegas_rad
    else:
        raise ValueError('Unknown mode_type {}: should be either \'axial\' or \'radial\''.format(mode_type))

    if len(used_ions) != 2: 
        raise ValueError('Two used ions in a chain must be specified')
    if any(i < 0 or i >= chain.n_ions for i in used_ions):
        raise ValueError("Ion indices out of range")

    eta1 = eta[used_ions[0], :]
    eta2 = eta[used_ions[1], :]

    return get_ppoly_pulse(eta1, eta2, omegas, *args, **kwargs)



def eff_omega_from_omega(x):
    return x*(j0(2*x) + jv(2, 2*x))


MAX_OMEGABYMU = 0.5818652242815963

class ExceedMaxOmegaError(ValueError):
    pass

def omega_from_eff_omega(y):
    if abs(y) > MAX_OMEGABYMU: # max value of eff_omega_from_omega
        raise ExceedMaxOmegaError('Equation for Omega from effective Omega cannot be solved')
    
    x = fsolve(lambda x: eff_omega_from_omega(x) - y, y)[0]
    return x


def Omega_max_abs(t_bp, Omega_mat):
    def spline_extrema(spline):
        return spline.derivative().roots()
    
    def spline_max_abs(spline):
        return np.max(np.abs(spline(spline_extrema(spline))))

    Omega = PPoly.construct_fast(Omega_mat, t_bp)
    return spline_max_abs(Omega)


def get_opt_Omega_mat(t_bp, mu, Omega_mat, N1):
    if Omega_max_abs(t_bp, Omega_mat)/mu > MAX_OMEGABYMU:
        raise ExceedMaxOmegaError('Omega(t)/mu exceeds MAX_OMEGABYMU')

    def divide_each_segment(t_bp, n_div):
        t_bp_1 = []
        for i in range(1, t_bp.size):
            t_bp_1.extend(np.linspace(t_bp[i-1], t_bp[i], n_div, endpoint=False))
        t_bp_1.append(t_bp[-1])
        return np.array(t_bp_1)

    t_bp_1 = divide_each_segment(t_bp, N1)
    
    Omegas = PPoly.construct_fast(Omega_mat, t_bp)
    Omega_eff_1_vals = np.array([Omegas(t) for t in t_bp_1])
    Omega_1_vals     = np.array([mu*omega_from_eff_omega(Omega/mu) for Omega in Omega_eff_1_vals])
    Omega_1_spline   = CubicSpline(t_bp_1, Omega_1_vals)
    Omega_1_mat      = Omega_1_spline.c

    return t_bp_1, Omega_1_mat


def get_nseg(shaping_type, n_df):
    if shaping_type in ['pointwise spline', '3']:
        return n_df + 1
    elif shaping_type in ['smoothed stepwise', '0s']:
        return 2*n_df + 1
    elif shaping_type in ['linear', '1']:
        return n_df + 1
    elif shaping_type in ['constant piecewise', '0']:
        return n_df
    else:
        raise ValueError('Unsupported shaping type')


def get_pmax(shaping_type):
    if shaping_type in ['pointwise spline', '3']:
        pmax = 3
    elif shaping_type in ['smoothed stepwise', '0s']:
        pmax = 3
    elif shaping_type in ['linear', '1']:
        pmax = 1
    elif shaping_type in ['constant piecewise', '0']:
        pmax = 0
    else:
        raise ValueError('Unsupported shaping type')

    return pmax
