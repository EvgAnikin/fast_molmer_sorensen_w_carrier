import math
import numpy as np

from scipy.interpolate import CubicSpline, PPoly

from . import pulse_shaping
from .ms_evo_funcs.ppoly_alpha_and_chi import *
from .ms_evo_funcs.num_ppoly_alpha_and_chi import num_alpha_and_chi_eqf, num_alpha_and_chi_eqf_w_carrier


class MSHandlerPPolyEqF:
    @classmethod 
    def get_from_chain(cls, chain, t_bp, Omega_mat, used_ions, mode_type, muby2pi, psibypi, anglebypi):
        mu = 2*math.pi*muby2pi
        psi = math.pi*psibypi
        angle = math.pi*anglebypi
        
        if mode_type == 'axial':
            eta = chain.eta_ax(angle)[used_ions, :]
            omegas = chain.omegas_ax
        elif mode_type == 'radial':
            eta = chain.eta_rad(angle)[used_ions, :]
            omegas = chain.omegas_rad
        else:
            raise ValueError('Unknown mode_type {}: should be either \'axial\' or \'radial\''.format(mode_type))

        return MSHandlerPPolyEqF(eta, omegas, mu, psi, t_bp, Omega_mat)

    @classmethod
    def get_from_chain_with_pulse_shaping(cls, chain, t_gate, thetabypi, n_df, shaping_type, smoothing_degree=0.1, carrier_opt=False, n_div=20, **attrs):
        t_bp, Omega_mat = pulse_shaping.get_ppoly_pulse_from_chain(chain, t_gate=t_gate, thetabypi=thetabypi, n_df=n_df, 
                                         shaping_type=shaping_type, smoothing_degree=smoothing_degree, **attrs)
        if carrier_opt:
            mu = attrs['muby2pi']*2*np.pi
            t_bp, Omega_mat = pulse_shaping.get_opt_Omega_mat(t_bp, mu, Omega_mat, n_div)

        return MSHandlerPPolyEqF.get_from_chain(chain, t_bp, Omega_mat, **attrs)

    @classmethod
    def get_with_pulse_shaping(cls, t_gate, eta, omegas, muby2pi, psibypi, thetabypi, n_df, shaping_type, smoothing_degree=0.1, 
            carrier_opt=False, n_div=20):
        mu  = 2*math.pi*muby2pi
        psi = math.pi*psibypi
        t_bp, Omega_mat = pulse_shaping.get_ppoly_pulse(eta[0,:], eta[1,:], omegas, muby2pi, t_gate, psibypi, thetabypi, n_df, shaping_type, smoothing_degree)

        if carrier_opt:
            t_bp, Omega_mat = pulse_shaping.get_opt_Omega_mat(t_bp, mu, Omega_mat, n_div)

        return MSHandlerPPolyEqF(eta, omegas, mu, psi, t_bp, Omega_mat)


    @classmethod
    def get_with_eff_omega(cls, handler, N_1):
        n_seg = len(handler.t_bp) - 1
        t_bp_1 = np.linspace(handler.t_i, handler.t_f, 1 + n_seg*N_1)
        mu = handler.mu
        psi = handler.psi
        Omega_eff_1_vals = np.array([handler.Omega(t) for t in t_bp_1])
        Omega_1_vals = np.array([handler.mu*pulse_shaping.eff_omega_from_omega(Omega/handler.mu) for Omega in Omega_eff_1_vals])
        Omega_1_spline = CubicSpline(t_bp_1, Omega_1_vals)
        Omega_1_mat = Omega_1_spline.c

        return  MSHandlerPPolyEqF(handler.eta, handler.omegas, mu, psi, t_bp_1, Omega_1_mat)


    @classmethod
    def get_from_eff_omega(cls, handler, N_1):
        n_seg = len(handler.t_bp) - 1
        t_bp_1 = np.linspace(handler.t_i, handler.t_f, 1 + n_seg*N_1)
        mu = handler.mu
        psi = handler.psi
        Omega_eff_1_vals = np.array([handler.Omega(t) for t in t_bp_1])
        Omega_1_vals = np.array([handler.mu*pulse_shaping.omega_from_eff_omega(Omega/handler.mu) for Omega in Omega_eff_1_vals])
        Omega_1_spline = CubicSpline(t_bp_1, Omega_1_vals)
        Omega_1_mat = Omega_1_spline.c

        return  MSHandlerPPolyEqF(handler.eta, handler.omegas, mu, psi, t_bp_1, Omega_1_mat)


    def __init__(self, eta, omegas, mu, psi, t_bp, Omega_mat):
        self.eta = eta
        self.omegas = omegas
        self.mu = mu
        self.psi = psi
        self.t_i  = t_bp[0]
        self.t_f  = t_bp[-1]

        self.n_used_ions = self.eta.shape[0]
        self.psis = psi*np.ones(self.n_used_ions)
        self.n_modes = self.eta.shape[1]
        self.n_seg = len(t_bp) - 1

        self.Omega_mat = Omega_mat
        self.t_bp = t_bp
        self.pmax = Omega_mat.shape[0] - 1

        if self.pmax > 3: 
            raise ValueError("zero dimension of Omega_mat exceeds 3: only polynomials up to 3rd order are implemented")
        if Omega_mat.shape[1] != self.n_seg: 
            raise ValueError("""the shape of Omega_mat should match the number of segments,
                                                     Omega_mat.shape: {} N_seg: {}""".format(Omega_mat.shape, self.N_seg))

        
        self.Omega = PPoly.construct_fast(Omega_mat, t_bp)
        self.create_carrier_phases_at_bp()
        self.create_alpha_at_bp()
        self.create_chi_at_bp()


    def create_carrier_phases_at_bp(self):
        self.carrier_phases_at_bp = ppoly_carrier_phase_eqf_at_bp(self.t_bp, self.Omega_mat, self.mu, self.psi)

    def create_alpha_at_bp(self):
        self.mode_alpha_at_bp = ppoly_alpha_eqf_at_bp(self.t_bp, self.Omega_mat, self.mu, self.psi, self.omegas)
        self.alpha_at_bp = self.eta[np.newaxis, :, :]*self.mode_alpha_at_bp[:, np.newaxis, :]

    def create_chi_at_bp(self):
        self.mode_chi_at_bp = ppoly_chi_eqf_at_bp(self.t_bp, self.Omega_mat, self.mu, self.psi, self.omegas, self.mode_alpha_at_bp)
        self.chi_at_bp = np.real(np.dot(self.eta[np.newaxis, :, :]*self.mode_chi_at_bp[:, np.newaxis, :], self.eta.T))

    def carrier_phase(self, t2, t1):
        return ppoly_carrier_phase_eqf(t2, t1, self.t_bp, self.Omega_mat, self.mu, self.psi, self.carrier_phases_at_bp)

    def alpha_func(self, t2, t1):
        mode_alpha = ppoly_alpha_eqf(t2, t1, self.t_bp, self.Omega_mat, self.mu, self.psi, self.omegas, self.mode_alpha_at_bp)
        return self.eta*mode_alpha[np.newaxis, :]

    def chi_func(self, t2, t1):
        mode_chi = ppoly_chi_eqf(t2, t1, self.t_bp, self.Omega_mat, self.mu, self.psi, self.omegas, self.mode_alpha_at_bp, self.mode_chi_at_bp)
        return (self.eta*mode_chi) @ self.eta.T

    def calculate_num_alpha_and_chi(self, n_vals):
        t_range = np.linspace(self.t_i, self.t_f, n_vals)

        alpha_and_chi = num_alpha_and_chi_eqf(n_vals, self.Omega_mat, self.t_bp, self.mu, self.psi, self.omegas)

        self.mode_alpha_num = alpha_and_chi[0]
        self.mode_chi_num = np.real(alpha_and_chi[1])

        self.t_range = t_range
        self.alpha_num = self.eta[np.newaxis, :, :]*self.mode_alpha_num[:, np.newaxis, :]
        self.chi_num = np.real(np.dot(self.eta[np.newaxis, :, :]*self.mode_chi_num[:, np.newaxis, :], self.eta.T))
        self.alpha_fin_num = self.alpha_num[-1, :, :]
        self.chi_fin_num = self.chi_num[-1, :, :]

    def calculate_num_alpha_and_chi_w_carrier(self, n_vals):
        t_range = np.linspace(self.t_i, self.t_f, n_vals)

        alpha_and_chi = num_alpha_and_chi_eqf_w_carrier(n_vals, self.Omega_mat, self.t_bp, self.mu, self.psi, self.omegas, self.carrier_phases_at_bp)

        self.mode_alpha_num_w_carrier = alpha_and_chi[0]
        self.mode_chi_num_w_carrier = np.real(alpha_and_chi[1])

        self.t_range = t_range
        self.alpha_num_w_carrier = self.eta[np.newaxis, :, :]*self.mode_alpha_num_w_carrier[:, np.newaxis, :]
        self.chi_num_w_carrier = np.real(np.dot(self.eta[np.newaxis, :, :]*self.mode_chi_num_w_carrier[:, np.newaxis, :], self.eta.T))
        self.alpha_fin_num_w_carrier = self.alpha_num_w_carrier[-1, :, :]
        self.chi_fin_num_w_carrier   = self.chi_num_w_carrier[-1, :, :]


    def to_dict(self):
        hdict = {}
    
        hdict['eta']     = self.eta.tolist()
        hdict['nus']  = (self.omegas/(2*math.pi)).tolist()
        hdict['muby2pi'] = self.mu/(2*math.pi)
        hdict['psis']     = self.psis.tolist()
        hdict['t_bp'] = self.t_bp.tolist()
        hdict['Omega_mat'] = np.real(self.Omega_mat).tolist()
    
        return hdict
