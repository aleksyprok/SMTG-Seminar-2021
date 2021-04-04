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

data_dir = 'u_perp_wave_kx=1_alpha=tan_1_4'
getgrid(DIR = data_dir + '/Data/')

file_number = len(glob.glob(data_dir + '/Data/ux1*'))

output_dir = data_dir + '/Figures/Line'
os.makedirs(output_dir, exist_ok = True)

n_vars = 5
var_list = ['ux1', 'bx1', 'b_par1', 'u_perp1', 'b_perp1']
title_list = [r'$u_x / u_0$', r'$b_x / b_0$', r'$b_{||} / b_0$', r'$u_\perp / u_0$', r'$b_\perp / b_0$']
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

y0 = 0
nz_plot = 256
plot_z = np.linspace(0, 4, nz_plot)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = 9.92 * 2
fig_size[1] = 3.59 * 2.5
fig.set_size_inches(fig_size)

for n in range(file_number):

    time = getstep(n, DIR = data_dir + '/Data/') * dt

    i = -1
    for var in var_list:

        i += 1

        coord1 = i % 3
        coord2 = i // 3
        ax = fig.add_subplot(2, 3, i + 1)

        get_xyz(var, DIR = data_dir + '/Data/')
        var1 = getdata(n, var, DIR = data_dir + '/Data/')

        var1_interp_func = interpolate.interp2d(y, z, var1, kind = 'cubic')
        plot_var1 = var1_interp_func(y0, plot_z)


        ax.plot(plot_z, plot_var1)
        ax.set_title(title_list[i])
        ax.set_ylim(min_var_array[i], max_var_array[i])
        ax.get_xaxis().set_visible(False)

        if i >= 2:
            ax.get_xaxis().set_visible(True)
            ax.set_xlabel(r'$z / L_0$')
        if i == 4:
            ax.text(1.5, 0.8, \
                	r'$t / t_0 = $' + '{:.2f}'.format(time), \
                    fontsize = 20, \
                	transform=ax.transAxes)
            ax.text(1.5, 0.1, \
                    r'$y = $' + '{:.1f}'.format(y0) + '\n' \
                	r'$k_x / k_{||} = $' + '{:.1f}'.format(kx / k_par) + '\n' + \
                    r'$\tan(\alpha) = $' + '{:.3f}'.format(np.tan(alpha)) + '\n' + \
                	r'$k_{||} = \pi / (2L_0)$' + '\n' \
                    r'$b_0 = B_0u_0 / v_{A0}$' + '\n' \
                    r'$t_0 = L_0 / v_{A0}$', \
                	transform=ax.transAxes)

        # if i == 1:
        #     ax.text(-0.4, 1.1, \
        #         	r'$k_x / k_{||} = $' + '{:.1f}'.format(kx / k_par) + '\n' + \
        #             r'$\tan(\alpha) = $' + '{:.3f}'.format(np.tan(alpha)), \
        #         	transform=ax.transAxes)
        # if i == 2:
        #     ax.text(0.3, 1.1, \
        #         	r'$t / t_0 = $' + '{:.2f}'.format(time), \
        #             fontsize = 20, \
        #         	transform=ax.transAxes)
        # if i == 3:
        #     ax.text(0.8, 1.1, \
        #         	r'$k_{||} = \pi / (2L_0)$' + '\n' \
        #             r'$b_0 = B_0u_0 / v_{A0}$', \
        #         	transform=ax.transAxes)

    fig.savefig(output_dir + '/' + '{:03d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()
