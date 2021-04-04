# coding: utf-8
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pyPLUTO as pp
from scipy import interpolate

frame=2258
L_scale=1e9
max_radius=14000
amr_level=3 #this seems to work well between resolution and time to plot the full domain

D = pp.pload(frame,datatype='hdf5',level=amr_level)
R, Theta=np.meshgrid(D.x1, D.x2)
R=R.T
Theta=Theta.T
idx=np.where( (R<max_radius))
X=R[idx]*np.sin(Theta[idx])
Y=R[idx]*np.cos(Theta[idx])
plot_x=np.linspace(X.min(), X.max(), num=1000)
plot_y=np.linspace(Y.min(), Y.max(), num=1000)

plot_X, plot_Y = np.meshgrid(plot_x, plot_y)
points=np.empty([X.size, 2])
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

del(R,Theta, plot_y, plot_x)

Z = interpolate.griddata(points, D.rho[idx].flatten(), (plot_X, plot_Y), method='nearest', rescale=True)

fig = plt.figure()
ax = plt.subplot(111)


im=ax.imshow(np.log10(Z), origin='lower', extent=[0,max_radius*L_scale, 0, max_radius*L_scale],cmap=plt.get_cmap('magma'))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.00)
cbar = plt.colorbar(im, cax=cax, ticklocation='right')
cbar.set_label(r'Log($\rho$)', size=10)

ax.ticklabel_format(style='scientific', useMathText=True)

ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
plt.savefig('16TI_pluto_%d.pdf'%(frame), bbox_inches='tight')
