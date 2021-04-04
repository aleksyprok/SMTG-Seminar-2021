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

data_dir = 'alfven_wave_data_and_figures'
getgrid(DIR = data_dir + '/Data/')

file_number = len(glob.glob(data_dir + '/Data/ux1*'))

output_dir = data_dir + '/Figures/Contour'
os.makedirs(output_dir, exist_ok = True)

var = 'ux1'

min_var1 = 0
max_var1 = 1
nlev = 50
levels = np.linspace(min_var1, max_var1, nlev)

ny_plot = 256
nz_plot = 256
plot_y = np.linspace(-5, 5, ny_plot)
plot_z = np.linspace(-5, 5, nz_plot)
plot_Y, plot_Z = np.meshgrid(plot_y, plot_z)

fig = plt.figure()

for n in range(file_number):

    ax = fig.add_subplot(111)

    get_xyz(var, DIR = data_dir + '/Data/')
    Y, Z = np.meshgrid(y, z)
    var1 = getdata(n, var, DIR = data_dir + '/Data/')

    var1_interp_func = interpolate.interp2d(y, z, var1, kind = 'cubic')
    plot_var1 = var1_interp_func(plot_y, plot_z)

    plot_var1 += (max_var1 - plot_var1) * (plot_var1 > max_var1)
    plot_var1 += (min_var1 - plot_var1) * (plot_var1 < min_var1)

    cs = ax.contourf(plot_Y, plot_Z, plot_var1, cmap = cm.jet, levels = levels)
    cs.set_clim([min_var1, max_var1])

    A = np.cos(alpha) * plot_Y - np.sin(alpha) * plot_Z
    fieldlines = plt.contour(plot_Y, plot_Z, A, colors = 'black')

    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    fig.savefig(output_dir + '/' + '{:03d}'.format(n) + '.png', bbox_inches='tight', dpi = 200)
    fig.clf()
