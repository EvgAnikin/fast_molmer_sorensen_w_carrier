import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import scipy.special as special


#In [18]: res = minimize(lambda x: -eff_omega(x), 1)
#
#In [19]: res.x
#Out[19]: array([0.92059315])
#
#In [20]: -res.fun
#Out[20]: 0.5818652242803048


MAX_OMEGAEFF = 0.5818652242803048
ARGMAX       = 0.92059315

def omega_eff_dm(x):
    return x*(special.jv(0, 2*x) + special.jv(2, 2*x))


if __name__ == '__main__':
    plt.figure(figsize=(5,2.9))
    matplotlib.rcParams.update({'font.size': 11})
    
    x = np.linspace(0, 3.5, 201)
    y = omega_eff_dm(x)

#    plt.plot([ARGMAX, ARGMAX], [-1, MAX_OMEGAEFF], linestyle='--', color='tab:blue')
    plt.plot([-1, ARGMAX], [MAX_OMEGAEFF, MAX_OMEGAEFF], linestyle='--', color='tab:orange')

    plt.text(0.06, 0.47, f'{MAX_OMEGAEFF:.4f}', color='tab:orange')
#    plt.text(1, -0.35, f'{ARGMAX:.4f}', color='tab:blue')

    plt.plot(x, y, linestyle='-', color='#AA1111', label=r'$x(J_0(2x) + J_2(2x))$')
    plt.ylim(-0.4, 0.7)
    plt.xlim(x[0], x[-1])
    plt.xlabel('$x = \Omega/\mu$')
    plt.ylabel('$\Omega_\mathrm{eff}/\mu$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
