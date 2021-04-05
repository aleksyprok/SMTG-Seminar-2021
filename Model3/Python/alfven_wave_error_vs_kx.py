import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size': 16})

def exp_cap(z):
	return np.exp(z * (z <= 100))

def ux_uniform(z):
	for n in range(3):
		if n == 0:
			ux  = ux0[n] * np.exp(1j * kz[n] * z)
		else:
			ux += ux0[n] * np.exp(1j * kz[n] * z)
	return ux

def ux_piecewise(z):
	for n in range(4):
		if n == 0:
			ux  = ux0_m[n] * exp_cap(1j * m_m[n] * z) * np.heaviside(-z, 0)
			ux += ux0_p[n] * exp_cap(1j * m_p[n] * z) * np.heaviside( z, 1)
		else:
			ux += ux0_m[n] * exp_cap(1j * m_m[n] * z) * np.heaviside(-z, 0)
			ux += ux0_p[n] * exp_cap(1j * m_p[n] * z) * np.heaviside( z, 1)
	return ux

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

L = 1
vA_p = 1.0
vA0 = vA_p
u0 = 1
omega = np.pi * vA_p / L

vA_m = 1e-2

n_kx = 1000
n_alpha = 3
n_vA = 3
kx_array = np.logspace(-1, 4, n_kx) * omega / vA_p
tan_alpha_array = np.array([1e-1, 1, 1e1])
alpha_array = np.arctan(tan_alpha_array)
vA_array = np.array([1e-1, 1e-2, 1e-3])
print(alpha_array / np.pi)

k_par_p = omega / vA_p
k_par_m = omega / vA_m

clrs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 1.5
fig.set_size_inches(fig_size)

ax = fig.add_subplot(121)

lines = [Line2D([0], [0], color = 'k', linestyle='--'), \
         Line2D([0], [0], color = 'k', linestyle=':' )]
labels = [r'$k_x$ = $k_{crit+}$', \
          r'$k_x$ = $k_{crit-}$']
lgd1 = ax.legend(lines, labels, bbox_to_anchor=(1.25,1.25))

sqrt_arr_m = np.zeros(n_alpha)
sqrt_arr_p = np.zeros(n_alpha)
for ialpha in range(n_alpha):

	alpha = alpha_array[ialpha]
	err_array = np.zeros(n_kx)
	ky = k_par_p * np.cos(alpha)

	for ikx in range(n_kx):

		kx = kx_array[ikx]

		[ux0, u_perp0, bx0, b_perp0, b_par0, kz] = get_uniform_coeffs(kx, ky, vA_m, alpha)

		[u_perp0_m, u_perp0_p, ux0_m, ux0_p, bx0_m, bx0_p, \
		 b_perp0_m, b_perp0_p, b_par0_m, b_par0_p, m_m, m_p] = get_piecewise_coeffs(kx, ky, vA_m, alpha)

		err_array[ikx] = np.abs((u_perp0[1] - u_perp0_p[1]) / u_perp0_p[0])

	sqrt_arr_p[ialpha] = np.real(np.sqrt(k_par_p ** 2 - ky ** 2 + 0j)) / k_par_p
	sqrt_arr_m[ialpha] = np.real(np.sqrt(k_par_m ** 2 - ky ** 2 + 0j)) / k_par_p
	ax.loglog(kx_array / k_par_p, err_array.real, color = clrs[ialpha])
	ax.set_ylim(1e-3, 1e-1)

lines = [Line2D([0], [0], color = 'tab:blue'), \
         Line2D([0], [0], color = 'tab:orange'), \
         Line2D([0], [0], color = 'tab:green')]
labels = [r'$\tan(\alpha)=10^{-1}$', r'$\tan(\alpha)=1$', r'$\tan(\alpha)=10^1$']
lgd = ax.legend(lines, labels, bbox_to_anchor=(0.3,-0.14))
ax.text(0.5,-0.325, \
    	r'$k_{||-} = 10^2 k_{||+}$', \
    	transform=ax.transAxes)
plt.gca().add_artist(lgd1)

yr = ax.get_ylim()
for ialpha in range(n_alpha):
    ln = ax.loglog([sqrt_arr_p[ialpha], sqrt_arr_p[ialpha]], [1e-12,1e12], linestyle = '--' , color = clrs[ialpha])
    ln = ax.loglog([sqrt_arr_m[ialpha], sqrt_arr_m[ialpha]], [1e-12,1e12], linestyle = ':' , color = 'k')
ax.set_ylim(yr)
ax.set_xlabel(r'$k_x / k_{||+}$')
ax.set_title(r'$\frac{\left|u_{\perp3} - u_{\perp3+}\right|}{\left|u_{\perp3+}\right|}$')

ax = fig.add_subplot(122)

alpha = np.pi / 4
ky = k_par_p * np.cos(alpha)
sqrt_arr_m = np.zeros(n_vA)
for ivA in range(n_vA):

	vA_m = vA_array[ivA]
	k_par_m = omega / vA_m
	err_array = np.zeros(n_kx)

	for ikx in range(n_kx):

		kx = kx_array[ikx]

		[ux0, u_perp0, bx0, b_perp0, b_par0, kz] = get_uniform_coeffs(kx, ky, vA_m, alpha)

		[u_perp0_m, u_perp0_p, ux0_m, ux0_p, bx0_m, bx0_p, \
		 b_perp0_m, b_perp0_p, b_par0_m, b_par0_p, m_m, m_p] = get_piecewise_coeffs(kx, ky, vA_m, alpha)

		err_array[ikx] = np.abs((u_perp0[1] - u_perp0_p[1]) / u_perp0_p[0])

	sqrt_arr_m[ivA] = np.real(np.sqrt(k_par_m ** 2 - ky ** 2 + 0j)) / k_par_p
	ln = ax.loglog(kx_array / k_par_p, err_array, color = clrs[ivA])
	clr = ln[0].get_color()
	yr = ax.get_ylim()

yr = ax.get_ylim()
ln = ax.loglog([np.real(np.sqrt(k_par_p ** 2 - ky ** 2 + 0j)) / k_par_p, \
	            np.real(np.sqrt(k_par_p ** 2 - ky ** 2 + 0j)) / k_par_p], \
	            [1e-12, 1e12], 'k--')
for ivA in range(n_vA):
	ln = ax.loglog([sqrt_arr_m[ivA], sqrt_arr_m[ivA]], [1e-12,1e12], linestyle = ':' , color = clrs[ivA])
ax.set_ylim(yr)
ax.set_xlabel(r'$k_x / k_{||+}$')
ax.set_title(r'$\frac{\left|u_{\perp3} - u_{\perp3+}\right|}{\left|u_{\perp3+}\right|}$')

lines = [Line2D([0], [0], color = 'tab:blue'), \
         Line2D([0], [0], color = 'tab:orange'), \
         Line2D([0], [0], color = 'tab:green')]
labels = [r'$k_{||-} = 10^{1}k_{||+}$', r'$k_{||-} = 10^{2}k_{||+}$', r'$k_{||-} = 10^{3}k_{||+}$']
lgd = ax.legend(lines, labels, bbox_to_anchor=(0.475,-0.14))
ax.text(0.5,-0.325, \
    	r'$\alpha = \pi / 2$', \
    	transform=ax.transAxes)

fig.savefig('Figures/alfven_wave_error_vs_kx.pdf', bbox_inches='tight')
