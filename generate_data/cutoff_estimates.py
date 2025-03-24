import itertools
import math
import numpy as np


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
