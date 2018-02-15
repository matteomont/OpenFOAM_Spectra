#!/usr/bin/env python


# Script that take a probe-surface file (x-z plane versus time)
# and generate premultiplied spectra, and saves spectra of velocity
# and pressure

#copyright by Matteo Montecchia

import os
import matplotlib.pyplot as plt 
import time as t
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab as p
import numpy as np
import numpy.matlib
from scipy.fftpack import fft,ifft
from scipy import signal
#from matplotlib import rc
#rc('text', usetex = True)


start_time = t.time()

velocity = input('Digit True if you want velocity spectra, False for pressure\n') 
filename = input('Digit file name to load\n')
n = np.loadtxt(filename)

print " --- file loading time = %s seconds ---" % (t.time() - start_time)



time = n[:,0]
Lz = 2*np.pi
Retau = 481.43642 
nu = 1.3764e-5
utau = Retau*nu

Nz = 128
if velocity:
	Ny = ((len(n[1,:])-1)/3)/Nz
else:
	Ny = (len(n[1,:])-1)/Nz

print Ny

yp = np.linspace(0,1,Ny)*Retau


if velocity:
	shift = 3
else:
	shift = 1
print shift
r=1
u = np.zeros((len(time),Ny,Nz))
for i in range(0,Ny):
	for j in range(0,Nz):
                u[:,i,j] = n[:,r]
                r= r+shift 

# -----mean velocity field------#
U = np.zeros((Ny,Nz)) 
for i in range(0,len(time)):
        U = U + u[i,:,:]
U = U/len(time)

# -----fluctuating field-------#
for i in range(0,Ny):
	for j in range(0,Nz):
		u[:,i,j] = u[:,i,j] - U[i,j]
# ---------- urms**2 --------------#
u2 = u**2
U2 = np.zeros((Ny,Nz))
for i in range(0,len(time)):
        U2 = U2 + u2[i,:,:]
U2 = U2/len(time)

Nzv = np.linspace(0,Nz,Nz)
kz=(2*np.pi/Lz)*Nzv
kzp=kz[0:Nz/2]/Retau

# ////////////////////////// --- Welch's method --- ///////////////////////////////////#
P = np.zeros((len(time),Ny,Nz/2+1))
for i in range(0,len(time)):
        fx,P[i,:,:] = signal.welch(u[i,:,:],scaling='spectrum')

P2 = np.zeros((Ny,Nz/2+1))
for i in range(0,Ny):
	for j in range(0,Nz/2+1):
        	P2[i,j] = np.average(P[:,i,j])

lambdazp = 2*np.pi/kz[0:Nz/2+1]*Retau

#premultiply P2
P2P = np.zeros((Ny,Nz/2+1))
for i in range(0,Ny):
	P2P[i,:] = kz[0:Nz/2+1]*P2[i,:]

#/////////////////////////////////////////////////////////////////////////////////////////#

#------------------- Premultiplied spectrum contour plot ---------------------------------#
plt.figure(0)
[lambdazpp,ypp] = np.meshgrid(lambdazp[0:Nz/2],yp)
#cp = plt.pcolormesh(lambdazp,yp,np.abs(P2P[:,:-1])/(utau**2), shading = 'gouraud')
if velocity:
	cp = plt.contourf(lambdazpp,ypp,np.abs(P2P[:,:-1])/(utau**2), shading = 'gouraud')
	cp2 = plt.contour(lambdazpp,ypp,np.abs(P2P[:,:-1])/(utau**2), shading = 'gouraud')
else:
	cp = plt.contourf(lambdazpp,ypp,np.abs(P2P[:,:-1])/(utau**4), shading = 'gouraud')
	cp2 = plt.contour(lambdazpp,ypp,np.abs(P2P[:,:-1])/(utau**4), shading = 'gouraud')

plt.colorbar()
plt.yscale('log')
plt.xscale('log')
plt.ylim([1,1e2])
plt.xlim([100,500])
plt.xlabel('$\\lambda_z^+$',fontsize=25)
plt.ylabel('$k_z E_{uu}^{1D} / u_\\tau^2$',fontsize=25)
plt.subplots_adjust(top = 0.95, bottom = 0.2, left = 0.2, right = 0.95)
plt.tick_params(axis='both',which='major',labelsize=20.0)
plt.grid()

#-----------------------------------------------------------------------------------------#

#------------------Surface 3D plot------------------------------------------#
fig = plt.figure(1)
ax = fig.gca(projection='3d')

KZP,YP = np.meshgrid(kzp,yp)
if velocity:
	surf = ax.plot_surface(KZP,YP,P2[:,:-1]/(utau**2),cmap = cm.coolwarm,
  	linewidth = 0,antialiased=False)
	ax.set_zlabel('$E_{uu}^+$',fontsize=20)
else:
	surf = ax.plot_surface(KZP,YP,P2[:,:-1]/(utau**4),cmap = cm.coolwarm,
  	linewidth = 0,antialiased=False)
	ax.set_zlabel('$E_{pp}^+$',fontsize=20)

fig.colorbar(surf,shrink = 0.5, aspect=5)
ax.set_xlabel('$k_z^+$',fontsize=20)
ax.set_ylabel('$y^+$',fontsize=20)

ax.xaxis.set_scale('log')
ax.yaxis.set_scale('log')
ax.zaxis.set_scale('log')
#---------------------------------------------------------------------------#
if velocity:
	# save spectrum
	np.savetxt('Euup.txt', P2[:,:-1]/(utau**2) , fmt=" %5.6f")
	# save premultiplied spectrum
	np.savetxt('kzEuup.txt', P2P[:,:-1]/(utau**2), fmt=" %5.6f")
	np.savetxt('kz.txt', kz, header = 'kz', fmt=" %5.6f")
	np.savetxt('yp.txt', yp, header = 'y+', fmt=" %5.6f")
else:
	# save spectrum
	np.savetxt('Eppp.txt', P2[:,:-1]/(utau**4) , fmt=" %5.6f")
	# save premultiplied spectrum
	#np.savetxt('kzppp.txt', P2P[:,:-1]/(utau**4), fmt=" %5.6f")
	

print " --- %s seconds ---" % (t.time() - start_time)
plt.show()


