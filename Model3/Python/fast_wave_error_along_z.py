import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

def exp_cap(z):
	return np.exp(z * (z <= 100))

def ux_uniform(z):
	for n in range(3):
		if n == 0:
			ux  = ux0[n] * np.exp(1j * kz[n] * z)
		else:
			ux += ux0[n] * np.exp(1j * kz[n] * z)
	return ux

def b_par_uniform(z):
	for n in range(3):
		if n == 0:
			b_par  = b_par0[n] * np.exp(1j * kz[n] * z) * (z >= 0)
		else:
			b_par += b_par0[n] * np.exp(1j * kz[n] * z) * (z >= 0)
	return b_par

def ux_piecewise(z):
	for n in range(4):
		if n == 0:
			ux  = ux0_m[n] * exp_cap(1j * m_m[n] * z) * np.heaviside(-z, 0)
			ux += ux0_p[n] * exp_cap(1j * m_p[n] * z) * np.heaviside( z, 1)
		else:
			ux += ux0_m[n] * exp_cap(1j * m_m[n] * z) * np.heaviside(-z, 0)
			ux += ux0_p[n] * exp_cap(1j * m_p[n] * z) * np.heaviside( z, 1)
	return ux

def b_par_piecewise(z):
	for n in range(4):
		if n == 0:
			b_par  = b_par0_m[n] * exp_cap(1j * m_m[n] * z) * (z <  0)
			b_par += b_par0_p[n] * exp_cap(1j * m_p[n] * z) * (z >= 0)
		else:
			b_par += b_par0_m[n] * exp_cap(1j * m_m[n] * z) * (z <  0)
			b_par += b_par0_p[n] * exp_cap(1j * m_p[n] * z) * (z >= 0)
	return b_par


def get_uniform_coeffs(kx, ky, vA_m, alpha):

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

	return [ux0, u_perp0, bx0, b_perp0, b_par0, kz]

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

output_dir = 'Figures/fast_wave_error_along_z'
os.makedirs(output_dir, exist_ok = True)

L0 = 1

vA_p = 1.0
vA0 = vA_p
vA_m = 1e-2 * vA_p
u0 = 1
omega = 2 * np.pi * vA_p / L0

k_par_p = omega / vA_p
k_par_m = omega / vA_m

alpha = np.pi / 6
kx = 1 * omega / vA0
ky = np.sin(alpha) * k_par_p

kx_crit_p = k_par_p * np.cos(alpha)
kx_crit_m = np.sqrt(k_par_m ** 2 - ky ** 2)

nz = 256
z_min = 0
z_max = L0
dz = (z_max - z_min) / (nz - 1)
z = np.linspace(z_min, z_max, nz)

nz_small = 256
z_small_min = -2 * np.pi / k_par_m
z_small_max =  2 * np.pi / k_par_m
dz_small = (z_small_max - z_small_min) / (nz_small - 1)
z_small = np.linspace(z_small_min, z_small_max, nz_small)

nz_rhs = 128
z_rhs_min = 0
z_rhs_max = 2 * np.pi / k_par_m
dz_rhs = (z_rhs_max - z_rhs_min) / (nz_rhs - 1)
z_rhs = np.linspace(z_rhs_min, z_rhs_max, nz_rhs)

n_kx = 101
kx_array = kx_crit_p * np.logspace(-1, 3, n_kx)

print(kx_array / kx_crit_p)
print(kx_array / kx_crit_m)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.2)

n = -1

for kx in kx_array:

    n += 1

    [ux0, u_perp0, bx0, b_perp0, b_par0, kz] = get_uniform_coeffs(kx, ky, vA_m, alpha)

    [u_perp0_m, u_perp0_p, ux0_m, ux0_p, bx0_m, bx0_p, \
     b_perp0_m, b_perp0_p, b_par0_m, b_par0_p, m_m, m_p] = get_piecewise_coeffs(kx, ky, vA_m, alpha)

    ax = fig.add_subplot(221)
    ax.plot(z_small / (2 * np.pi / k_par_m), b_par_piecewise(z_small).real)
    ax.plot(z_rhs / (2 * np.pi / k_par_m), b_par_uniform(z_rhs).real)
    ax.set_ylim(-1, 1)

    ax = fig.add_subplot(222)
    ax.plot(z, b_par_piecewise(z).real)
    ax.plot(z, b_par_uniform(z).real)
    ax.set_ylim(-1, 1)

    ax = fig.add_subplot(223)
    ax.plot(z_small / (2 * np.pi / k_par_m), b_par_piecewise(z_small).imag)
    ax.plot(z_rhs / (2 * np.pi / k_par_m), b_par_uniform(z_rhs).imag)
    ax.set_ylim(-1, 1)

    ax = fig.add_subplot(224)
    ax.plot(z, b_par_piecewise(z).imag)
    ax.plot(z, b_par_uniform(z).imag)
    ax.set_ylim(-1, 1)

    fig.savefig(output_dir + '/' + '{:04d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()


# plt.show()
