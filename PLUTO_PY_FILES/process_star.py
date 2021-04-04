import matplotlib.pyplot as plt;

plt.ion()
import numpy as np
from scipy.interpolate import interp1d
import os
m_p=1.6726231e-24 #grams
c_light=2.99792458e10 #cm/s

#define sca;es for normalization here
L_scale=1e9 #cm
dens_scale=1 #gr/cm^3
pres_scale=dens_scale*c_light**2 #g/cm/s^2

#define the coordinate system to be used in PLUTO
coord='Spherical' #spherical for r and theta, assume symmetry about origin so ignore phi
#coord='Cartesian'

r_solar = 6.96e10  # cm
r_16ti = 4.07e10  # cm
r_16oi=4.45936e10 #r_solar*0.7
r_35ob=3.92e10 #r_solar*0.6

#define the stellar profile to use
star='16TI'
r_star=r_16ti


# read in the data
data = np.genfromtxt(star+'.presn', skip_header=2) #data from https://2sn.org/GRB2/

# find the middle radius of each shell
avg_r = 0.5 * (data[:, 2][:-1] + data[:, 2][1:])

# find the mass in each shell
shell_mass = np.diff(data[:, 1])

# find the volume of each shell
vol = (4 / 3) * np.pi * (-data[:, 2][:-1] ** 3 + data[:, 2][1:] ** 3)

# calculate the density in each shell and get the other relevant parameters
dens = shell_mass / vol

#get the pressure in each cell
pres=0.5*(data[:-1, 6]+data[1:, 6])

def expand1dGrid(avg_r, dens, coord):

    # interpolate the log(density) and normalized radius data
    f = interp1d(avg_r / r_solar, np.log10(dens), kind='cubic')

    # argument to f is radius/r_solar and the output is log10(density)
    # The interpolation for 16TI is only god until r_16ti/r_solar=0.58477, where r_16ti=4.07e10 cm

    # can check the interpolation by plotting
    plt.figure()
    plt.plot(avg_r / r_solar, np.log10(dens), label='data')
    x = np.linspace(avg_r.min(), 1 * 6.96e10)
    plt.plot(x / r_solar, f(x / r_solar), label='fit')
    plt.legend()

    # create grid and fill it with the values
    if coord == 'Cartesian':
        x = np.linspace(avg_r.min(), 1.2 * r_star, num=1000)
        y = x.copy()
    if coord == 'Spherical':
        x = np.linspace(avg_r.min(), 1.2 * r_star, num=1000)
        y = np.linspace(0, np.pi/2, num=1000)

    X, Y = np.meshgrid(x, y)

    # save values b/c that may be useful later
    X_old = X.copy()
    Y_old = Y.copy()

    # make it cell centered values
    X = X + np.diff(x)[0] / 2
    Y = Y + np.diff(y)[0] / 2

    # save the values at the cell  centers that are within the stars' radius
    dens_profile = np.zeros_like(X)
    if coord == 'Cartesian':
        idx = np.where(np.sqrt(X ** 2 + Y ** 2) < r_star)
        dens_profile[idx] = f(np.sqrt(X[idx] ** 2 + Y[idx] ** 2) / r_solar)
    else:
        idx = np.where(X  < r_star)
        dens_profile[idx] = f(X[idx] / r_solar)


    #plt.figure()
    #plt.contourf(X,Y,dens_profile)
    #plt.colorbar(pad=0)

    # the values outside of r_16ti are wrong so we need to do a spline in this range and then refill these grid points
    # with the correct values
    idx = np.where(avg_r >= r_star)[0]  # this may just be the very last value in the array which is the only node outside r_16ti,
    # but need 2 points there so take r_16ti and that same last point to do the interpolation
    f2 = np.polyfit(avg_r[idx[0] - 1:] / r_solar, np.log10(dens[idx[0] - 1:]), 1)
    f_f_2 = np.poly1d(f2)

    # modify grid
    if coord == 'Cartesian':
        idx = np.where(np.sqrt(X ** 2 + Y ** 2) >= r_star)
        dens_profile[idx] = f_f_2(np.sqrt(X[idx] ** 2 + Y[idx] ** 2) / r_solar)
    else:
        idx = np.where(X  >= r_star)
        dens_profile[idx] = f_f_2(X[idx]  / r_solar)


    # replot to check it
    plt.figure()
    if coord=='Cartesian':
        x=X
        y=Y
    else:
        x=X*np.cos(Y)
        y=X*np.sin(Y)
    plt.contourf(X,Y,dens_profile)#-2*np.log10(c_light))
    plt.colorbar(pad=0)
    #plt.plot(X,Y, marker='.', color='k', linestyle='none')
    #plt.plot(X_old, Y_old, marker='.', color='r', linestyle='none')

    return dens_profile, X, Y, X_old, Y_old

dens_profile, X, Y, X_old, Y_old=expand1dGrid(avg_r, dens, coord)
pres_profile, X, Y, X_old, Y_old=expand1dGrid(avg_r, pres, coord)

#normalize the lengths and the density and the pressure by the same units that pluto uses
if coord=='Cartesian':
    X, Y, X_old, Y_old=X/L_scale, Y/L_scale, X_old/L_scale, Y_old/L_scale
else:
    X, Y, X_old, Y_old = X / L_scale, Y , X_old / L_scale, Y_old
#dens_profile=dens_profile/dens_scale
#pres_profile=pres_profile/pres_scale

# save the reformatted data
if not os.path.exists('STAR_FORMATTED'):
    os.makedirs('STAR_FORMATTED')

if  not os.path.exists('STAR_FORMATTED/'+star+'_star_grid.out'):

    X_L=X_old[0,:]
    X_R=X_old[0,:]+np.diff(X_old[0,:])[0]
    i=np.arange(1,np.size(X_L)+1)

    Y_L=Y_old[:,0]
    Y_R=Y_old[:,0]+np.diff(Y_old[:,0])[0]
    j=np.arange(1,np.size(Y_L)+1)

    if coord=='Cartesian':
        header='# GEOMETRY:   CARTESIAN\n'
    else:
        header='# GEOMETRY:   SPHERICAL\n'

    f=open('STAR_FORMATTED/'+star+'_star_grid.out','a+')
    f.write('%s'%(header))
    f.write('%d\n'%(np.size(X_L)))
    np.savetxt(f, np.transpose((i,X_L, X_R)), fmt="%d\t%e\t%e")

    f.write('%d\n'%(np.size(Y_L)))
    np.savetxt(f, np.transpose((j,Y_L, Y_R)), fmt="%d\t%e\t%e")

    #do third dimension because pluto expects it
    f.write('%d\n'%(1))
    f.write('%s'%("1  0.0  1.0\n"))

    f.close()

    #now save the normalized variables that we want pluto to read in to interpolate
    ((10**dens_profile)/dens_scale).flatten().tofile('STAR_FORMATTED/'+star+'_star_dens.dbl')
    ((10**pres_profile)/pres_scale).flatten().tofile('STAR_FORMATTED/'+star+'_star_pres.dbl')

else:
    print('There is already a processed file in the directory. Please delete it and then rerun the python program')

