import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D
import os
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams.update({'font.size': 16})

def reflected_alfven(y, z, t, var):
    if var == 'ux':
        return ux0[1] * np.exp(1j * (ky * y + kz[1] * z + omega * t))
    elif var == 'u_perp':
        return u_perp0[1] * np.exp(1j * (ky * y + kz[1] * z + omega * t))
    elif var == 'bx':
        return bx0[1] * np.exp(1j * (ky * y + kz[1] * z + omega * t))
    elif var == 'b_perp':
        return b_perp0[1] * np.exp(1j * (ky * y + kz[1] * z + omega * t))
    elif var == 'b_par':
        return b_par0[1] * np.exp(1j * (ky * y + kz[1] * z + omega * t))
    else:
        print('Invalid string.')

L0 = 1

ny = 256
y_min = -1.5 * L0
y_max =  1.5 * L0
dy = (y_max - y_min) / (ny - 1)
y = np.linspace(y_min, y_max, ny)

nz = 256
z_min = 0
z_max = 2 * L0
dz = (z_max - z_min) / (nz - 1)
z = np.linspace(z_min, z_max, nz)

Y, Z = np.meshgrid(y, z)

vA0 = 1
alpha = np.pi / 6
omega = 2 * np.pi
k_par = omega / vA0
ky = np.sin(alpha) * k_par
kx_crit = np.cos(alpha) * k_par
u0 = 1

output_dir = 'Figures'
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

n_kx = 1001
kx_array = k_par * np.logspace(-2, 2, n_kx)

n_alpha = 3
alpha_array = np.array([np.arctan(1e1), np.arctan(1e0), np.arctan(1e-1)])

ux1_array = np.zeros(n_kx, dtype = complex)
u_perp1_array = np.zeros(n_kx, dtype = complex)
u_amp1_array = np.zeros(n_kx)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

for alpha in alpha_array:

    ky = np.sin(alpha) * k_par
    kx_crit = np.cos(alpha) * k_par
    n = -1
    print(kx_crit / k_par)

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

        ux1_array[n] = ux0[1]
        u_perp1_array[n] = u_perp0[1]
        u_amp1_array[n] = np.sqrt(np.abs(ux0[1]) ** 2 + np.abs(u_perp0[1]) ** 2)

    ax1.plot(kx_array / k_par, u_amp1_array)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$k_x / k_{||}$')
    ax1.set_title(r'$\langle|\mathbf{u}|\rangle / u_0$')

    # ax.plot(kx_array / kx_crit, ux1_array.real)
    # ax.plot(kx_array / kx_crit, ux1_array.imag)
    # ax.plot(kx_array / kx_crit, u_perp1_array.real)
    # ax.plot(kx_array / kx_crit, u_perp1_array.imag)
    ax2.plot(kx_array / k_par, np.arctan2(u_perp1_array.real, ux1_array.real) / np.pi)
    ax2.set_ylim(-0.5, 0)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'$k_x / k_{||}$')
    ax2.set_title(r'$\phi / \pi$')

line1 = Line2D([1],[1], color = "tab:blue", \
               label = r'$\tan(\alpha) = $' + '{:.1f}'.format(np.tan(alpha_array[0])))
line2 = Line2D([1],[1], color = "tab:orange", \
               label = r'$\tan(\alpha) = $' + '{:.2f}'.format(np.tan(alpha_array[1])))
line3 = Line2D([1],[1], color = "tab:green", \
               label = r'$\tan(\alpha) = $' + '{:.2f}'.format(np.tan(alpha_array[2])))
ax1.add_line(line1)
ax1.add_line(line2)
ax1.add_line(line3)
handles, labels = ax1.get_legend_handles_labels()
lgd = ax1.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.075,1.2))

fig.savefig(output_dir + '/' + 'reflect_alfven_amp_polar.pdf', bbox_inches='tight')
plt.show()
