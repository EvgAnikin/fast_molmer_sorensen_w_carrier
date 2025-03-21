import collections
import datetime
import itertools
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import timeit
import sys

from scipy.interpolate import CubicSpline, PPoly


def constant_piecewise_matrix(t_bp, pmax=0):
    n_bp = len(t_bp) 
    T = np.zeros((pmax+1, n_bp - 1, n_bp - 1))

    for k in range(1, n_bp):
        dt = t_bp[k] - t_bp[k-1]
        T[pmax, k-1, k-1] = 1/dt

    return T

def linear_interp_matrix(t_bp, pmax=1):
    n_bp = len(t_bp) 
    T = np.zeros((pmax+1, n_bp - 1, n_bp - 2))

    for k in range(1, n_bp - 1):
        dt_l = t_bp[k] - t_bp[k-1]
        dt_r = t_bp[k+1] - t_bp[k]

        T[pmax-1, k-1, k-1] = 1/dt_l
        T[pmax, k-1, k-1] = 0
        T[pmax-1,   k, k-1] = -1/dt_r
        T[pmax,   k, k-1] = 1

    return T


def pointwise_spline_matrix(t_bp):
    n_bp = len(t_bp) 
    T = np.zeros((4, n_bp - 1, n_bp - 2))

    for k in range(1, n_bp - 1):
        values = np.zeros(n_bp)
        values[k] = 1
        basis_spline = CubicSpline(t_bp, values, bc_type='clamped')
        T[:, :, k - 1] = basis_spline.c

    return T


def breakpoints_for_smoothed_piecewise(n_plat, t_plat, t_trans, t0=0):
    t_bp = [t0]
    for i in range(n_plat):
        t_bp.append(t_bp[-1] + t_trans)
        t_bp.append(t_bp[-1] + t_plat)

    t_bp.append(t_bp[-1] + t_trans)
    return np.array(t_bp)


def smoothed_piecewise_matrix(t_bp):
    n_seg = len(t_bp) - 1
    n_plateaus = int((n_seg - 1)/2)

    T = np.zeros((4, n_seg, n_plateaus))
    for k in range(n_plateaus):
        # Transition function between segments: 3*x**2 - 2*x**3, where x = (t - t[k])/(t[k+1] - t[k])
        dt_l = t_bp[2*k+1] - t_bp[2*k]
        T[0, 2*k, k] = -2/dt_l**3
        T[1, 2*k, k] = 3/dt_l**2

        T[3, 2*k + 1, k] = 1


        dt_r = t_bp[2*k+3] - t_bp[2*k+2]
        T[0, 2*k+2, k] = 2/dt_r**3
        T[1, 2*k+2, k] = -3/dt_r**2
        T[3, 2*k+2, k] = 1

    return T


def ppoly_norm_matrix(t_bp, pmax):
    n_seg = len(t_bp) - 1
    N = np.zeros((pmax+1, n_seg, pmax+1, n_seg))

    def segment_norm_matrix(dt):
        p_range = np.arange(0, pmax+1)
        p_sum = (p_range[:, np.newaxis] + p_range[np.newaxis, :]) 
        return 1/(p_sum+1)*dt**(p_sum + 1)

    for i in range(n_seg):
        N[:, i, :, i] = segment_norm_matrix(t_bp[i+1] - t_bp[i])

    return N


def ppoly_norm_matrix_reversed(t_bp, pmax):
    n_seg = len(t_bp) - 1
    N = np.zeros((pmax+1, n_seg, pmax+1, n_seg))

    def segment_norm_matrix(dt):
#        p_range = np.arange(0, pmax+1)
        p_range = np.arange(pmax, -1, -1)
        p_sum = (p_range[:, np.newaxis] + p_range[np.newaxis, :]) 
        return 1/(p_sum+1)*dt**(p_sum + 1)

    for i in range(n_seg):
        N[:, i, :, i] = segment_norm_matrix(t_bp[i+1] - t_bp[i])

    return N



def pulse_square_integral(t_bp, Omega_mat):
    pmax = Omega_mat.shape[0] - 1
    norm_matrix = ppoly_norm_matrix(t_bp, pmax)
    Omega_mat_rv = Omega_mat[::-1, ::-1]
    return np.tensordot(Omega_mat_rv, np.tensordot(norm_matrix, Omega_mat_rv))
