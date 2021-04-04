import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile
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

# var = 'u_perp1'
var = input('Enter variable name: ')
file_number = len(glob.glob('Data/' + var + '*'))

output_dir = 'Figures/Surface/' + var + '/'
os.makedirs(output_dir, exist_ok = True)

# for n in range(file_number):
#
# 	if n % 100 == 0: print(n)
#
# 	time = getstep(n) * dt
# 	var1 = getdata(n, var)
# 	if n == 0:
# 		max_var1 = np.max(var1)
# 		min_var1 = np.max(var1)
# 	else:
# 		max_var1 = np.max([max_var1, np.max(var1)])
# 		min_var1 = np.min([min_var1, np.min(var1)])

min_var1 = -1
max_var1 = 1

print("max_var1 = " + str(max_var1))
print("min_var1 = " + str(min_var1))

fig = plt.figure()

for n in range(file_number):

	time = getstep(n) * dt

	ax = fig.add_subplot(111, projection = '3d')

	get_xyz(var)
	Y, Z = np.meshgrid(y, z)
	var1 = getdata(n, var)
	surf = ax.plot_surface(Y, Z, var1, cmap = cm.cool)
	ax.set_zlim([min_var1, max_var1])
	surf.set_clim([min_var1, max_var1])
	ax.set_xlabel('y')
	ax.set_ylabel('z')
	ax.set_title(var + ', t = ' + '{:.2f}'.format(time))

	fig.savefig(output_dir + '/' + '{:03d}'.format(n) + '.png', bbox_inches='tight')
	fig.clf()
