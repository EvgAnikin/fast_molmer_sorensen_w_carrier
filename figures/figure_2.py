import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, LinearSegmentedColormap
import math
from scipy.interpolate import PPoly
from timslib.ion_crystals.ion_chain import IonChain
from timslib.ion_crystals.normal_modes import nu_ax_from_nu_rad_min
from timslib.ms_pulse_shaping.ms_handler import MSHandlerPPolyEqF as MSHandler

def plot_pictures(chain, handler):
    n_carrier = 2501
    matrix = chain.eta_rad(shaping_attrs["anglebypi"]*np.pi)  
    freqs = chain.omegas_rad/(2*np.pi*1e6)
    t_range = np.linspace(0, shaping_attrs["t_gate"], n_carrier)

    chis_opt = []
    chis_no_opt = []

    n = 2
    aph_opt = np.zeros((n,2,len(t_range)))
    aph_no_opt = np.zeros((n,2,len(t_range)))
    
    shaping_attrs["carrier_opt"] = True
    chain = IonChain(**chain_attrs)
    handler = MSHandler.get_from_chain_with_pulse_shaping(chain, **shaping_attrs)
    handler.calculate_num_alpha_and_chi_w_carrier(n_carrier)

    Om_opt = handler.Omega(t_range)

    for j in range(n):
        aph_opt[j,0,:] = handler.alpha_num_w_carrier[:,0,j+3].real
        aph_opt[j,1,:] = handler.alpha_num_w_carrier[:,0,j+3].imag
        chis_opt = handler.chi_num_w_carrier[:,0,1]
    
    shaping_attrs["carrier_opt"] = False
    handler = MSHandler.get_from_chain_with_pulse_shaping(chain, **shaping_attrs)
    handler.calculate_num_alpha_and_chi_w_carrier(n_carrier)

    Om_no_opt = handler.Omega(t_range)

    for j in range(n):
        aph_no_opt[j,0,:] = handler.alpha_num_w_carrier[:,0,j+3].real
        aph_no_opt[j,1,:] = handler.alpha_num_w_carrier[:,0,j+3].imag
        chis_no_opt = handler.chi_num_w_carrier[:,0,1]
    
    def plot_frequency_lines(ax, freqs):
        for freq in freqs:
            ax.hlines(freq, xmin=0, xmax=0.9, color='blue')
        ax.get_xaxis().set_visible(False)
        ax.set_ylim([freqs.min() - 0.05, freqs.max() + 0.05])
        ax.set_yticks(freqs)
        ax.set_yticklabels([f'{freq:.2f}' for freq in freqs])
        ax.set_ylabel('$\omega_m/(2\pi)$ (MHz)')
        ax.grid(True, axis='y')
        ax.text(-0.15, 1.1, 'b)', transform=ax.transAxes, fontsize=12, va='top')
#        ax.yaxis.tick_right()
#        ax.yaxis.set_label_position("right")
    
    def plot_matrix(ax, matrix):
        cmap = plt.get_cmap('seismic')
        x_range = np.arange(0.5, 5.5 + 1e-10)
        y_range = np.arange(0.5, 5.5 + 1e-10)
        cax = ax.pcolormesh(x_range, y_range, matrix.T, cmap=cmap, vmin=-0.08, vmax=0.08)

        rect = matplotlib.patches.Rectangle((1.5,0.53), 2, 4.92, color='#60ff46', 
                fill=False, linewidth=3)
        ax.add_patch(rect)

        cbar = fig.colorbar(cax, ax=ax, location='right', shrink=0.8)
        cbar.set_label(r'$\tilde\eta_{im}^\mathrm{rad}$')
        ax.set_aspect('equal')
        ax.set_xlabel('Ion index in the chain')
        ax.set_ylabel('Mode index')
        ax.set_xticks(range(1,6))
        ax.set_yticks(range(1,6))
#        ax_twin = ax.twiny()
#        ax_twin.set_xticks([2,3])
        ax.set_xlim(0.5,5.5)
#        ax_twin.set_xlim(0.5,5.5)
        ax.set_ylim(0.5, 5.5)
        ax.text(-0.15, 1.1, 'c)', transform=ax.transAxes, fontsize=12, va='top')
        ax.text(2, 5.7, '1', ha='center', fontsize=12)
        ax.text(3, 5.7, '2', ha='center', fontsize=12)
    
    def plot_omega(ax, t_range, Om_opt, Om_no_opt):
        ax.plot(t_range*1e6, Om_opt/1e3/(2*np.pi),    label=r'$\Omega_\mathrm{tr}(t)$')
        ax.plot(t_range*1e6, Om_no_opt/1e3/(2*np.pi), linestyle='--', label=r'$\Omega_\mathrm{lin}(t)$')
        ax.set_ylabel(r'$\Omega(t)/(2\pi)$, kHz')
        ax.set_xlabel(r't, $\mu s$')
        ax.legend(loc='lower right')
        ax.text(-0.15, 1.1, 'd)', transform=ax.transAxes, fontsize=12, va='top')
        ax.margins(0)
        ax.set_ylim(0, max(Om_opt)/1e3*1.1/(2*np.pi))
    
    def plot_alpha_lines(ax, aph_opt, aph_no_opt, n):
        for j in range(n):
            ax.plot(aph_opt[j][0], aph_opt[j][1], label=fr'$\alpha_{{1{j+4}}}$ tr')

        for j in range(n):
            ax.plot(aph_no_opt[j][0], aph_no_opt[j][1], label=fr'$\alpha_{{1{j+4}}}$ lin', linestyle='--')


    def plot_alpha(ax, aph_opt, aph_no_opt, n):
        plot_alpha_lines(ax, aph_opt, aph_no_opt, n)
        axins = ax.inset_axes([0.21, 0.68, 0.35, 0.35], xlim=(-0.02,  0.02), ylim=(-0.02, 0.02))#, xticklabels=[], yticklabels=[])
        plot_alpha_lines(axins, aph_opt, aph_no_opt, n)
        axins.xaxis.set_ticks_position('top')

        ax.set_ylabel(r'$Im(\alpha)$')
        ax.set_xlabel(r'$Re(\alpha)$')
        ax.text(-0.15, 1.1, 'e)', transform=axs[1,1].transAxes, fontsize=12, va='top')
        ax.legend(loc='lower right')

        ax.set_yticks([-0.2, 0, 0.2])
        ax.grid()

        ax.indicate_inset_zoom(axins, edgecolor="black")

        alpha_noopt_fin_0 = aph_no_opt[0][0][-1] + 1j*aph_no_opt[0][1][-1]
        alpha_noopt_fin_1 = aph_no_opt[1][0][-1] + 1j*aph_no_opt[1][1][-1]
        axins.scatter([alpha_noopt_fin_0.real], [alpha_noopt_fin_0.imag], color='tab:green', marker='D', zorder=3, s=40)
        axins.scatter([alpha_noopt_fin_1.real], [alpha_noopt_fin_1.imag], color='tab:red', marker='D', zorder=3, s=40)

        axins.grid()

    
    def plot_chi(ax):
        ax.plot(t_range*1e6,chis_opt/np.pi, label=r'tr')
        ax.plot(t_range*1e6, chis_no_opt/np.pi, label=r'lin', linestyle='--')
        #ax.set_title(r'$\chi_{01}(t)/\pi$', fontsize=7)
        ax.set_xlabel(r't, $\mu s$')
        ax.set_ylabel(r'$\chi_{12}(t)/\pi$')
        ax.legend(loc='lower right')
        ax.text(-0.15, 1.1, 'f)', transform=ax.transAxes, fontsize=12, va='top')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.margins(0)
        ax.set_ylim(ymax=max(chis_opt/np.pi)*1.1)

    def plot_eq_pos(ax):
        eq_pos = chain.eq_pos
        x_range = np.linspace(eq_pos[0] - 1e-5, eq_pos[-1] + 1e-5, 501)
        y_range = np.linspace(-1e-5, 1e-5, 101)
        x, y = np.meshgrid(x_range, y_range)
        r = 1.3e-6
        total_intensity = sum(np.exp(-np.sqrt(((ion_x - x) ** 2 + y ** 2)) / r) for ion_x in eq_pos)
        ax.set_xlabel(r'$\mu m$')
        ax.pcolormesh(x_range * 1e6, y_range * 1e6, total_intensity, cmap=plt.get_cmap('Blues').reversed(), vmax=0.9)
        ax.text(-0.15, 1.1, 'a)', transform=ax.transAxes, fontsize=12, va='top')
        ax.text(eq_pos[1]*1e6, 4, '1', fontsize=14, color='white', ha='center')
        ax.text(eq_pos[2]*1e6, 4, '2', fontsize=14, color='white', ha='center')
        ax.set_aspect('equal')
        ax.set_yticks([])

    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    plot_eq_pos(axs[0, 0])
    plot_frequency_lines(axs[0, 1], freqs)
    plot_matrix(axs[0, 2], matrix)
    plot_omega( axs[1, 0], t_range, Om_opt, Om_no_opt)
    plot_alpha(axs[1,1], aph_opt, aph_no_opt, n)
    plot_chi(   axs[1, 2])
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 11})

    n_ions = 5
    nu_rad = 1e6
    nu_rad_min = 0.75e6
    nu_ax = nu_ax_from_nu_rad_min(nu_rad_min, nu_rad, n_ions)

    chain_attrs =  {
        "ion_type": "Ca40",
        "nu_ax": nu_ax,
        "nu_rad": nu_rad,
        'n_ions': n_ions
    }

    n1 = 1
    n2 = 2

    shaping_attrs = {
        "anglebypi": 0.5,
        "mode_type": "radial",
        "t_gate":  0.00004174150891,
        "muby2pi":  1033765.2642760169 ,
        "psibypi": 0.5,
        "thetabypi": 0.25,
        "shaping_type": "3",
        "used_ions": [n1, n2],
        "n_df": 2 * n_ions + 1,
        "carrier_opt": True,
    }
    chain = IonChain(**chain_attrs)
    handler = MSHandler.get_from_chain_with_pulse_shaping(chain, **shaping_attrs)
    plot_pictures(chain, handler)
