import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile
from scipy import interpolate
import glob
import os

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
           alpha, rho_bcc, rho_bbb, rho_ccb, rho_cbc

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

file_number = len(glob.glob('Data/' + var + '*'))

output_dir = 'Figures/Line/' + var
os.makedirs(output_dir, exist_ok = True)

min_var1 = -1
max_var1 = 1
nlev = 50
levels = np.linspace(min_var1, max_var1, nlev)

nz_plot = 256
y0 = 0
plot_z = np.linspace(0, 4, nz_plot)

fig = plt.figure()
# fig_size = fig.get_size_inches()
# fig_size[0] = 9.92 * 2.5
# fig_size[1] = 3.59 * 2.5
# fig.set_size_inches(fig_size)

for n in range(file_number):

    time = getstep(n) * dt

    ax = fig.add_subplot(111)

    time = getstep(n) * dt

    get_xyz(var)

    var1 = getdata(n, var)

    var1_interp_func = interpolate.interp2d(y, z, var1, kind = 'cubic')
    plot_var1 = var1_interp_func(0, plot_z)

    ax.plot(plot_z, plot_var1, '+-')
    ax.set_xlabel('z')
    ax.set_title(var + ', t = ' + '{:.2f}'.format(time))
    ax.set_ylim(min_var1, max_var1)

    fig.savefig(output_dir + '/' + '{:03d}'.format(n) + '.png', bbox_inches='tight')
    fig.clf()
