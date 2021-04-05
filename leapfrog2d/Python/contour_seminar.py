import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import os
plt.rcParams.update({'font.size': 16})
plt.rcParams['contour.negative_linestyle'] = 'solid'

def getdata(filenumber, variable, DIR = ' '):

    if DIR == ' ':
        filename = 'Data/' + variable + '{:03d}'.format(filenumber) + '.dat'
    else:
        filename = DIR + variable + '{:03d}'.format(filenumber) + '.dat'

    file = FortranFile(filename, 'r')
    ny = file.read_ints(np.int32)[0]
    nz = file.read_ints(np.int32)[0]
    return file.read_reals(float).reshape((nz,ny), order="C")

def getgrid(DIR = ' '):

    global t_max, dt, nx, ny, nz, xb, yb, zb, xc, yc, zc, \
           alpha, k_par, kx

    if DIR == ' ':
    	filename = 'Data/grid.dat'
    else:
    	filename = DIR + 'grid.dat'

    file = FortranFile(filename, 'r')

    t_max = file.read_reals(float)[0]
    dt = file.read_reals(float)[0]
    ny = file.read_ints(np.int32)[0]
    nz = file.read_ints(np.int32)[0]
    yb = file.read_reals(float)
    zb = file.read_reals(float)
    yc = file.read_reals(float)
    zc = file.read_reals(float)
    alpha = file.read_reals(float)[0]
    k_par = file.read_reals(float)[0]
    kx = file.read_reals(float)[0]

def get_xyz(variable, DIR = ' '):

	global x, y, z

	if DIR == ' ':
		getgrid()
	else:
		getgrid(DIR)

	if variable == 'ux1':
		y = yc
		z = zc
	elif variable == 'ux2':
		y = yb
		z = zb
	elif variable == 'u_perp1':
		y = yc
		z = zb
	elif variable == 'u_perp2':
		y = yb
		z = zc
	elif variable == 'bx1':
		y = yc
		z = zb
	elif variable == 'bx2':
		y = yb
		z = zc[1:-1]
	elif variable == 'b_perp1':
		y = yc
		z = zc[1:-1]
	elif variable == 'b_perp2':
		y = yb
		z = zb
	elif variable == 'b_par1':
		y = yc
		z = zc[1:-1]
	elif variable == 'b_par2':
		y = yb
		z = zb
	else:
		print('Error')

def getstep(filenumber, DIR = ' '):
	if DIR == ' ':
		filename = 'Data/step' '{:03d}'.format(filenumber) + '.dat'
	else:
		filename = DIR + 'step' + '{:03d}'.format(filenumber) + '.dat'
	file = FortranFile(filename, 'r')
	return file.read_ints(np.int32)[0]

def get_max_min_val(var):

    for n in range(file_number):

    	if n % 10 == 0: print(n)

    	var1 = getdata(n, var)
    	if n == 0:
    		max_var1 = np.max(var1)
    		min_var1 = np.max(var1)
    	else:
    		max_var1 = np.max([max_var1, np.max(var1)])
    		min_var1 = np.min([min_var1, np.min(var1)])

    if abs(max_var1) < 1e-3: max_var1 =  1e-3
    if abs(min_var1) < 1e-3: min_var1 = -1e-3

    return [max_var1, min_var1]

data_dir = 'ux_wave_kx=0_alpha=tan_1_4'
getgrid(DIR = data_dir + '/Data/')

file_number = len(glob.glob(data_dir + '/Data/ux1*'))

output_dir = data_dir + '/Figures/Contour'
os.makedirs(output_dir, exist_ok = True)

n_vars = 5
var_list = ['ux1', 'u_perp1', 'bx1', 'b_perp1', 'b_par1']
title_list = [r'$u_x / u_0$', r'$u_\perp / u_0$', r'$b_x / b_0$', r'$b_\perp / b_0$', r'$b_{||} / b_0$']
max_var_array = np.zeros(n_vars, dtype = float)
min_var_array = np.zeros(n_vars, dtype = float)

n = -1
for var in var_list:
    n += 1
    [max_var1, min_var1] = get_max_min_val(var)
    max_var_array[n] = max_var1
    min_var_array[n] = min_var1
max_var_array = np.ones(n_vars)
min_var_array = -1 * np.ones(n_vars)

nlev = 25
levels_array = np.zeros((n_vars, nlev))
for n in range(n_vars):
    levels_array[n,:] += np.linspace(min_var_array[n], max_var_array[n], nlev)
#
ny_plot = 128
nz_plot = 256
plot_y = np.linspace(-1, 1, ny_plot)
plot_z = np.linspace(0, 4, nz_plot)
plot_Y, plot_Z = np.meshgrid(plot_y, plot_z)
#
fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)
# plt.subplots_adjust(left=0.125, bottom=0.1, right=0.5, top=0.9, wspace=0.2, hspace=0.2)

for n in range(file_number):

    time = getstep(n, DIR = data_dir + '/Data/') * dt

    i = -1
    for var in var_list:

        i += 1

        ax = fig.add_subplot(1, n_vars, i + 1)

        get_xyz(var, DIR = data_dir + '/Data/')
        Y, Z = np.meshgrid(y, z)
        var1 = getdata(n, var, DIR = data_dir + '/Data/')

        var1_interp_func = interpolate.interp2d(y, z, var1, kind = 'cubic')
        plot_var1 = var1_interp_func(plot_y, plot_z)

        plot_var1 += (max_var_array[i] - plot_var1) * (plot_var1 > max_var_array[i])
        plot_var1 += (min_var_array[i] - plot_var1) * (plot_var1 < min_var_array[i])

        cs = ax.contourf(plot_Y, plot_Z, plot_var1, cmap = cm.jet, levels = levels_array[i,:])
        cs.set_clim([min_var_array[i], max_var_array[i]])

        A = np.cos(alpha) * plot_Y - np.sin(alpha) * plot_Z
        fieldlines = plt.contour(plot_Y, plot_Z, A, colors = 'black')

        ax.set_xlabel(r'$y / L_0$')
        ax.set_ylabel(r'$z / L_0$')
        ax.set_title(title_list[i])
        ax.get_yaxis().set_visible(False)
        ax.set_aspect('equal')
        if i == 0:
            ax.get_yaxis().set_visible(True)
        if i == 1:
            ax.text(-0.4, 1.12, \
                	r'$k_x / k_{||} = $' + '{:.1f}'.format(kx / k_par) + '\n' + \
                    r'$\tan(\alpha) = $' + '{:.3f}'.format(np.tan(alpha)), \
                	transform=ax.transAxes)
        if i == 2:
            ax.text(0.3, 1.12, \
                	r'$t / t_0 = $' + '{:.2f}'.format(time), \
                    fontsize = 20, \
                	transform=ax.transAxes)
        if i == 3:
            ax.text(0.8, 1.1, \
                	r'$k_{||} = \pi / (2L_0)$' + '\n' \
                    r'$b_0 = B_0u_0 / v_{A0}$' + '\n' \
                    r'$t_0 = L_0 / v_{A0}$', \
                	transform=ax.transAxes)
        if i == n_vars - 1:
            # cb_ax = fig.add_axes([.91,.124,.01,.754])
            cb_ax = fig.add_axes([.91,.19,.01,.6])
            fig.colorbar(cs, orientation='vertical',cax=cb_ax)

    fig.savefig(output_dir + '/' + '{:03d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()

# n = 20
#
# time = getstep(n, DIR = data_dir + '/Data/') * dt
#
# i = -1
# for var in var_list:
#
#     i += 1
#
#     ax = fig.add_subplot(1, n_vars, i + 1)
#
#     get_xyz(var, DIR = data_dir + '/Data/')
#     Y, Z = np.meshgrid(y, z)
#     var1 = getdata(n, var, DIR = data_dir + '/Data/')
#
#     var1_interp_func = interpolate.interp2d(y, z, var1, kind = 'cubic')
#     plot_var1 = var1_interp_func(plot_y, plot_z)
#
#     plot_var1 += (max_var_array[i] - plot_var1) * (plot_var1 > max_var_array[i])
#     plot_var1 += (min_var_array[i] - plot_var1) * (plot_var1 < min_var_array[i])
#
#     cs = ax.contourf(plot_Y, plot_Z, plot_var1, cmap = cm.jet, levels = levels_array[i,:])
#     cs.set_clim([min_var_array[i], max_var_array[i]])
#
#     A = np.cos(alpha) * plot_Y - np.sin(alpha) * plot_Z
#     fieldlines = plt.contour(plot_Y, plot_Z, A, colors = 'black')
#
#     ax.set_xlabel(r'$y / L_0$')
#     ax.set_ylabel(r'$z / L_0$')
#     ax.set_title(title_list[i])
#     ax.get_yaxis().set_visible(False)
#     ax.set_aspect('equal')
#     if i == 0: ax.get_yaxis().set_visible(True)
#     if i == n_vars - 1:
#         # cb_ax = fig.add_axes([.91,.124,.01,.754])
#         cb_ax = fig.add_axes([.91,.19,.01,.6])
#         fig.colorbar(cs, orientation='vertical',cax=cb_ax)
#
#
# plt.show()
