import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams.update({'font.size': 16})

def reflected_fast(y, z, var):
    if var == 'ux':
        return ux0[2] * np.exp(1j * (ky * y + kz[2] * z))
    elif var == 'u_perp':
        return u_perp0[2] * np.exp(1j * (ky * y + kz[2] * z))
    elif var == 'bx':
        return bx0[2] * np.exp(1j * (ky * y + kz[2] * z))
    elif var == 'b_perp':
        return b_perp0[2] * np.exp(1j * (ky * y + kz[2] * z))
    elif var == 'b_par':
        return b_par0[2] * np.exp(1j * (ky * y + kz[2] * z))
    else:
        print('Invalid string.')

L0 = 1

nz = 256
z_min = 0
z_max = 2 * L0
dz = (z_max - z_min) / (nz - 1)
z = np.linspace(z_min, z_max, nz)

vA0 = 1
alpha = np.pi / 6
omega = 2 * np.pi
k_par = omega / vA0
ky = np.sin(alpha) * k_par
kx_crit = k_par * np.cos(alpha)
# kx = 1e-5
# kx = 0.1 * k_par
u0 = 1

output_dir = 'Figures/fast_wave_along_z'
os.makedirs(output_dir, exist_ok = True)

nt = 101
T = 2 * np.pi / omega
t_min = 0
t_max = 10 * T
t_array = np.linspace(t_min, t_max, nt)

nlev = 25
min_val = -1
max_val = 1
levels = np.linspace(min_val, max_val, nlev)

n_vars = 5
var_list = ['ux', 'bx', 'b_par', 'u_perp', 'b_perp']
title_list = [r'$u_x / u_0$', r'$b_x / b_0$', r'$b_{||} / b_0$', r'$u_\perp / u_0$', r'$b_\perp / b_0$']

kx_min = 0.1 * kx_crit
kx_max = 10 * kx_crit
# kx_array = np.arange(kx_min, kx_max + kx_min / 2, kx_min)
# n_kx = len(kx_array)
n_kx = 101
kx_array = kx_crit * np.logspace(-1, 1, n_kx)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)

n = -1

for kx in kx_array:

    n += 1

    kz0 = omega / vA0 / np.cos(alpha)
    kz = np.array([
    				 kz0 - ky * np.tan(alpha), \
    				-kz0 - ky * np.tan(alpha), \
    				 1j * np.sqrt(kx ** 2 + ky ** 2 - omega ** 2 / vA0 ** 2 + 0j), \
    				])

    nabla_par0  = 1j * (ky * np.sin(alpha) + kz * np.cos(alpha))
    nabla_perp0 = 1j * (ky * np.cos(alpha) - kz * np.sin(alpha))
    L0 = nabla_par0 ** 2 + omega ** 2 / vA0 ** 2
    ux_hat = -1j * kx * nabla_perp0 / (L0 - kx ** 2)

    u_perp0 = np.zeros(3, dtype=complex)
    u_perp0[0] =  u0
    u_perp0[1] = -u0 * (ux_hat[0] - ux_hat[2]) / (ux_hat[1] - ux_hat[2])
    u_perp0[2] =  u0 * (ux_hat[0] - ux_hat[1]) / (ux_hat[1] - ux_hat[2])

    ux0 = ux_hat * u_perp0
    bx0 = nabla_par0 * ux0 / (1j * omega)
    b_perp0 = nabla_par0 * u_perp0 / (1j * omega)
    b_par0 = -(1j * kx * ux0 + nabla_perp0 * u_perp0) / (1j * omega)

    i = -1
    for var in var_list:

        i += 1

        ax = fig.add_subplot(2, 3, i + 1)

        ax.plot(z, reflected_fast(0, z, var).real, label = 'Real part')
        ax.plot(z, reflected_fast(0, z, var).imag, label = 'Imag part')

        ax.set_title(title_list[i])

        ax.get_xaxis().set_visible(False)

        ax.set_ylim(min_val, max_val)

        if i >= 2:
            ax.get_xaxis().set_visible(True)
            ax.set_xlabel(r'$z / L_0$')

        if i == n_vars - 1:
            ax.text(1.5, 0.8, \
                	r'$k_x / k_{crit} = $' + '{:.2f}'.format(kx / kx_crit), \
                    fontsize = 20, \
                	transform=ax.transAxes)
            ax.text(1.25, 0.4, \
                	r'$k_{crit} = k_{||}\cos(\alpha)$'  '\n' + \
                    r'$\alpha = $' + '{:.3f}'.format(alpha / np.pi) + r'$\pi$', \
                	transform=ax.transAxes)
            ax.text(1.75, 0.4, \
                	r'$k_{||} = 2\pi / L_0$' + '\n' \
                    r'$b_0 = B_0u_0 / v_{A0}$', \
                	transform=ax.transAxes)
            handles, labels = ax.get_legend_handles_labels()
            lgd = ax.legend(handles, labels, bbox_to_anchor=(1.9,0.25))


    fig.savefig(output_dir + '/' + '{:04d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()

# plt.show()
