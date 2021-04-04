import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams.update({'font.size': 16})

def combined_wave(y, z, t, var):
    for n in range(4):
        if var == 'ux':
            if n == 0:
                wave  = ux0_p[0] * np.exp(1j * (ky * y + m_p[0] * z + omega * t)) * (z >= 0) + \
                        ux0_m[0] * np.exp(1j * (ky * y + m_m[0] * z + omega * t)) * (z  < 0)
            else:
                wave += ux0_p[n] * np.exp(1j * (ky * y + m_p[n] * z + omega * t)) * (z >= 0) + \
                        ux0_m[n] * np.exp(1j * (ky * y + m_m[n] * z + omega * t)) * (z  < 0)
        elif var == 'u_perp':
            if n == 0:
                wave  = u_perp0_p[0] * np.exp(1j * (ky * y + m_p[0] * z + omega * t)) * (z >= 0) + \
                        u_perp0_m[0] * np.exp(1j * (ky * y + m_m[0] * z + omega * t)) * (z  < 0)
            else:
                wave += u_perp0_p[n] * np.exp(1j * (ky * y + m_p[n] * z + omega * t)) * (z >= 0) + \
                        u_perp0_m[n] * np.exp(1j * (ky * y + m_m[n] * z + omega * t)) * (z  < 0)
        elif var == 'bx':
            if n == 0:
                wave  = bx0_p[0] * np.exp(1j * (ky * y + m_p[0] * z + omega * t)) * (z >= 0) + \
                        bx0_m[0] * np.exp(1j * (ky * y + m_m[0] * z + omega * t)) * (z  < 0)
            else:
                wave += bx0_p[n] * np.exp(1j * (ky * y + m_p[n] * z + omega * t)) * (z >= 0) + \
                        bx0_m[n] * np.exp(1j * (ky * y + m_m[n] * z + omega * t)) * (z  < 0)
        elif var == 'b_perp':
            if n == 0:
                wave  = b_perp0_p[0] * np.exp(1j * (ky * y + m_p[0] * z + omega * t)) * (z >= 0) + \
                        b_perp0_m[0] * np.exp(1j * (ky * y + m_m[0] * z + omega * t)) * (z  < 0)
            else:
                wave += b_perp0_p[n] * np.exp(1j * (ky * y + m_p[n] * z + omega * t)) * (z >= 0) + \
                        b_perp0_m[n] * np.exp(1j * (ky * y + m_m[n] * z + omega * t)) * (z  < 0)
        elif var == 'b_par':
            if n == 0:
                wave  = b_par0_p[0] * np.exp(1j * (ky * y + m_p[0] * z + omega * t)) * (z >= 0) + \
                        b_par0_m[0] * np.exp(1j * (ky * y + m_m[0] * z + omega * t)) * (z  < 0)
            else:
                wave += b_par0_p[n] * np.exp(1j * (ky * y + m_p[n] * z + omega * t)) * (z >= 0) + \
                        b_par0_m[n] * np.exp(1j * (ky * y + m_m[n] * z + omega * t)) * (z  < 0)
        else:
            print('Invalid string.')
    return wave

def get_piecewise_coeffs(kx, ky, vA_m, alpha):

	kz_m = omega / vA_m / np.cos(alpha)
	kz_p = omega / vA_p / np.cos(alpha)

	m_m = np.array([
					 kz_m - ky * np.tan(alpha), \
					-kz_m - ky * np.tan(alpha), \
					 1j * np.sqrt(kx ** 2 + ky ** 2 - omega ** 2 / vA_m ** 2 + 0j), \
					-1j * np.sqrt(kx ** 2 + ky ** 2 - omega ** 2 / vA_m ** 2 + 0j), \
					])
	m_p = np.array([
					 kz_p - ky * np.tan(alpha), \
					-kz_p - ky * np.tan(alpha), \
					 1j * np.sqrt(kx ** 2 + ky ** 2 - omega ** 2 / vA_p ** 2 + 0j), \
					-1j * np.sqrt(kx ** 2 + ky ** 2 - omega ** 2 / vA_p ** 2 + 0j), \
					])

	nabla_par_m  = 1j * (ky * np.sin(alpha) + m_m * np.cos(alpha))
	nabla_par_p  = 1j * (ky * np.sin(alpha) + m_p * np.cos(alpha))
	nabla_perp_m = 1j * (ky * np.cos(alpha) - m_m * np.sin(alpha))
	nabla_perp_p = 1j * (ky * np.cos(alpha) - m_p * np.sin(alpha))

	L_m = nabla_par_m ** 2 + omega ** 2 / vA_m ** 2
	L_p = nabla_par_p ** 2 + omega ** 2 / vA_p ** 2

	ux_hat_m = -1j * kx * nabla_perp_m / (L_m - kx ** 2)
	ux_hat_p = -1j * kx * nabla_perp_p / (L_p - kx ** 2)

	# Coefficent matrix
	aa = np.array([
					[         ux_hat_m[0],          ux_hat_m[3],          -ux_hat_p[1],          -ux_hat_p[2]], \
					[m_m[0] * ux_hat_m[0], m_m[3] * ux_hat_m[3], -m_p[1] * ux_hat_p[1], -m_p[2] * ux_hat_p[2]], \
					[                   1,                    1,                    -1,                    -1], \
					[              m_m[0],               m_m[3],               -m_p[1],               -m_p[2]]
					])
	bb = u0 * np.array([
						         ux_hat_p[0], \
						m_p[0] * ux_hat_p[0], \
						                   1, \
						              m_p[0], \
						])
	xx = np.linalg.solve(aa, bb)

	u_perp0_m = np.zeros(4, dtype=complex)
	u_perp0_p = np.zeros(4, dtype=complex)
	u_perp0_m[0] = xx[0]
	u_perp0_m[3] = xx[1]
	u_perp0_p[0] = u0
	u_perp0_p[1] = xx[2]
	u_perp0_p[2] = xx[3]


	ux0_m = ux_hat_m * u_perp0_m
	ux0_p = ux_hat_p * u_perp0_p

	bx0_m = nabla_par_m * ux0_m / (1j * omega)
	bx0_p = nabla_par_p * ux0_p / (1j * omega)

	b_perp0_m = nabla_par_m * u_perp0_m / (1j * omega)
	b_perp0_p = nabla_par_p * u_perp0_p / (1j * omega)

	b_par0_m = -(1j * kx * ux0_m + nabla_perp_m * u_perp0_m) / (1j * omega)
	b_par0_p = -(1j * kx * ux0_p + nabla_perp_p * u_perp0_p) / (1j * omega)

	return [u_perp0_m, u_perp0_p, ux0_m, ux0_p, bx0_m, bx0_p, \
	        b_perp0_m, b_perp0_p, b_par0_m, b_par0_p, m_m, m_p]

L0 = 1

ny = 256
y_min = -1.5 * L0
y_max =  1.5 * L0
dy = (y_max - y_min) / (ny - 1)
y = np.linspace(y_min, y_max, ny)

nz = 256
z_min = -L0
z_max =  L0
dz = (z_max - z_min) / (nz - 1)
z = np.linspace(z_min, z_max, nz)

Y, Z = np.meshgrid(y, z)

vA_p = 1.0
vA0 = vA_p
vA_m = 0.25 * vA_p
u0 = 1
omega = 2 * np.pi * vA_p / L0

k_par_p = omega / vA_p
k_par_m = omega / vA_m

alpha = np.pi / 6
kx = 1 * omega / vA0
ky = np.sin(alpha) * k_par_p


[u_perp0_m, u_perp0_p, ux0_m, ux0_p, bx0_m, bx0_p, \
 b_perp0_m, b_perp0_p, b_par0_m, b_par0_p, m_m, m_p] = get_piecewise_coeffs(kx, ky, vA_m, alpha)

output_dir = 'Figures/combined_wave'
os.makedirs(output_dir, exist_ok = True)

nt = 101
T = 2 * np.pi / omega
t_min = 0
t_max = 10 * T
t_array = np.linspace(t_min, t_max, nt)

nlev = 25
min_val = -1.4
max_val =  1.4
levels = np.linspace(min_val, max_val, nlev)

n_vars = 5
var_list = ['ux', 'bx', 'b_par', 'u_perp', 'b_perp']
title_list = [r'$u_x / u_0$', r'$b_x / b_0$', r'$b_{||} / b_0$', r'$u_\perp / u_0$', r'$b_\perp / b_0$']

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)

print(b_perp0_m[0])
print(b_perp0_m[1])

for n in range(nt):

    t = t_array[n]

    i = -1
    for var in var_list:

        i += 1

        ax = fig.add_subplot(2, 3, i + 1)

        cs = ax.contourf(Y, Z,
                         combined_wave(Y, Z, t, var).real, \
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
            ax.text(1.25, 0.35, \
                	r'$v_{A-} / v_{A+} = $' + '{:.2f}'.format(vA_m / vA_p) + '\n' + \
                	r'$k_x / k_{||+} = $' + '{:.1f}'.format(kx / k_par_p) + '\n' + \
                    r'$\alpha = $' + '{:.3f}'.format(alpha / np.pi) + r'$\pi$', \
                	transform=ax.transAxes)
            ax.text(1.75, 0.35, \
                	r'$k_{||+} = 2\pi / L_0$' + '\n' \
                    r'$b_0 = B_0u_0 / v_{A+}$' + '\n' \
                    r'$T = 2\pi / \omega$', \
                	transform=ax.transAxes)

        if i == n_vars - 1:
            cb_ax = fig.add_axes([.66,.15,.25,.015])
            cbar = fig.colorbar(cs, orientation = 'horizontal', cax=cb_ax)
            cbar.ax.tick_params(labelsize=12)


    fig.savefig(output_dir + '/' + '{:04d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()
