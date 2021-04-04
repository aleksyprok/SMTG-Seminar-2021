import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
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
kx = 1 * omega / vA0
ky = np.sin(alpha) * k_par
u0 = 1

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

output_dir = 'Figures/reflected_alfven'
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

print(np.amax(np.abs([ux0[1], u_perp0[1], bx0[1], b_perp0[1], b_par0[1]])))

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)

for n in range(nt):

    t = t_array[n]

    i = -1
    for var in var_list:

        i += 1

        ax = fig.add_subplot(2, 3, i + 1)

        cs = ax.contourf(Y, Z,
                         reflected_alfven(Y, Z, t, var).real, \
                         cmap = cm.jet, \
                         levels = levels, \
                         )
        cs.set_clim([min_val, max_val])
        ax.set_aspect('equal')
        ax.set_title(title_list[i])

        A = np.cos(alpha) * Y - np.sin(alpha) * Z
        fieldlines = plt.contour(Y, Z, A, colors = 'black')

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        if i >= 2:
            ax.get_xaxis().set_visible(True)
            ax.set_xlabel(r'$y / L_0$')

        if i % 3 == 0:
            ax.get_yaxis().set_visible(True)
            ax.set_ylabel(r'$z / L_0$')

        if i == n_vars - 1:
            ax.text(1.5, 0.8, \
                	r'$t / T = $' + '{:.2f}'.format(t / T), \
                    fontsize = 20, \
                	transform=ax.transAxes)
            ax.text(1.3, 0.325, \
                	r'$k_x / k_{||} = $' + '{:.1f}'.format(kx / k_par) + '\n' + \
                    r'$\alpha = $' + '{:.3f}'.format(alpha / np.pi) + r'$\pi$' '\n', \
                	transform=ax.transAxes)
            ax.text(1.7, 0.35, \
                	r'$k_{||} = 2\pi / L_0$' + '\n' \
                    r'$b_0 = B_0u_0 / v_{A0}$' + '\n' \
                    r'$T = 2\pi / \omega$', \
                	transform=ax.transAxes)

        if i == n_vars - 1:
            cb_ax = fig.add_axes([.66,.15,.25,.015])
            cbar = fig.colorbar(cs, orientation = 'horizontal', cax=cb_ax)
            cbar.ax.tick_params(labelsize=12)


    fig.savefig(output_dir + '/' + '{:04d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()
