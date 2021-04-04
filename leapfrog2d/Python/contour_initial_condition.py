import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile
from scipy import interpolate
import glob
import os
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

getgrid()

var = 'u_perp1'
data_dir = 'u_perp_wave_kx=10_alpha=atan_1_4'

file_number = len(glob.glob(data_dir + '/Data/ux1*'))

output_dir = 'Figures/'
os.makedirs(output_dir, exist_ok = True)

min_var1 = -1
max_var1 = 1
nlev = 25
levels = np.linspace(min_var1, max_var1, nlev)

ny_plot = 512
nz_plot = 256
plot_y = np.linspace(yb[0], yb[-1], ny_plot)
plot_z = np.linspace(zb[0], zb[-1], nz_plot)
plot_Y, plot_Z = np.meshgrid(plot_y, plot_z)

fig = plt.figure()

n = 0

ax = fig.add_subplot(111)

time = getstep(n, DIR = data_dir + '/Data/') * dt

get_xyz(var)
Y, Z = np.meshgrid(y, z)
var1 = getdata(n, var, DIR = data_dir + '/Data/')

var1_interp_func = interpolate.interp2d(y, z, var1, kind = 'cubic')
plot_var1 = var1_interp_func(plot_y, plot_z)

plot_var1 += (max_var1 - plot_var1) * (plot_var1 > max_var1)
plot_var1 += (min_var1 - plot_var1) * (plot_var1 < min_var1)


cs = ax.contourf(plot_Y, plot_Z, plot_var1, cmap = cm.jet, levels = levels)

A = np.cos(alpha) * plot_Y - np.sin(alpha) * plot_Z
fieldlines = plt.contour(plot_Y, plot_Z, A, colors = 'black')

ax.plot([-1, -1], [0, 4], 'b', linewidth = 3)
ax.plot([-1, 1], [4, 4], 'b', linewidth = 3)
ax.plot([1, 1], [0, 4], 'b', linewidth = 3)

cs.set_clim([min_var1, max_var1])
cbar = plt.colorbar(cs)
ax.set_xlabel(r'$y / L_0$')
ax.set_ylabel(r'$z / L_0$')
ax.set_title(r'$u_\perp / u_0$', fontsize = 16)
ax.set_aspect('equal')

ax.text(0.05, 1.05, \
        r'$t / t_0 = $' + '{:.2f}'.format(time) + '\n' + \
    	r'$k_x / k_{||} = $' + '{:.1f}'.format(kx / k_par) + '\n' + \
        r'$\tan(\alpha) = $' + '{:.3f}'.format(np.tan(alpha)), \
    	transform=ax.transAxes)
ax.text(0.75, 1.05, \
    	r'$k_{||} = \pi / (2L_0)$' + '\n' \
        r'$b_0 = B_0u_0 / v_{A0}$' + '\n' \
        r'$t_0 = L_0 / v_{A0}$', \
    	transform=ax.transAxes)

fig.savefig(output_dir + 'initial_condition.png', bbox_inches='tight', dpi = 200)
fig.clf()
