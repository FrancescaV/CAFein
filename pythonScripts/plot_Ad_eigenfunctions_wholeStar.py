import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
from math import sqrt
from math import sin
from math import pow
from math import fabs
matplotlib.rcParams.update({'font.size': 10})

#################################
## reading the constants
#################################
constants = "../constants.dat"

constantsVec   = np.loadtxt(constants, skiprows=0, usecols=(1,))

pi            = constantsVec[0]
G             = constantsVec[1]
cm_per_Rsun   = constantsVec[2]
g_per_Msun    = constantsVec[3]
g_per_Mearth  = constantsVec[4]
g_per_Mj      = constantsVec[5]
sec_per_Year  = constantsVec[6]
cm_per_Rj     = constantsVec[7]
cm_per_Rearth = constantsVec[8]
cm_per_AU     = constantsVec[9]
ergPerSec_per_Lsun = constantsVec[10]
BoltzmannK    = constantsVec[11]
AMU           = constantsVec[12]
sigma	      = constantsVec[13]
c_light_cgs   = constantsVec[14]

sec_per_day = 60.0*60.0*24.0
#################################
xmin = 0.
xmax = 1.1

basicAddress1 = "../output/"
addressFiles1 = basicAddress1+"/eigenfunctions.dat"
plotsAddress = basicAddress1+"plots/"

basicAddress2 = "../input/test_suite/polytrope_n3_eigenfunctions/Output/"
addressFiles2 = basicAddress2+"/eigenfunctions.dat"
#################################
# plot properties
#################################
#################################
fig_width_cm = 21                         # A4 page
fig_height_cm = 29.7
inches_per_cm = 1.0/2.54              # Convert cm to inches
fig_width = fig_width_cm * inches_per_cm # width in inches
fig_height = fig_height_cm * inches_per_cm       # height in inches
fig_size = [fig_width, fig_height]
##################################
##################################
PSfile = plotsAddress+"Eigenfunctions_"+str(xmin)+"_to_"+str(xmax)+".pdf"
plt.close('all')
fig, ((ax1yr), (ax2yr), (ax3yr), (ax4yr)) = plt.subplots(nrows=4, ncols=1, sharex = "all")
fig.set_size_inches(fig_size)
plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, wspace=0.3, hspace=0.2)


addressFiles = addressFiles1
rRel = np.loadtxt(addressFiles, skiprows=2, usecols=(0,))
y1r  = np.loadtxt(addressFiles, skiprows=2, usecols=(2,))
y2r  = np.loadtxt(addressFiles, skiprows=2, usecols=(3,))
y3r  = np.loadtxt(addressFiles, skiprows=2, usecols=(4,))
y4r  = np.loadtxt(addressFiles, skiprows=2, usecols=(5,))

data_1 = np.column_stack((rRel,y1r, y2r, y3r, y4r))
data_1 = data_1[np.argsort(data_1[:, 0])]

addressFiles = addressFiles2
rRel = np.loadtxt(addressFiles, skiprows=2, usecols=(0,))
y1r  = np.loadtxt(addressFiles, skiprows=2, usecols=(2,))
y2r  = np.loadtxt(addressFiles, skiprows=2, usecols=(3,))
y3r  = np.loadtxt(addressFiles, skiprows=2, usecols=(4,))
y4r  = np.loadtxt(addressFiles, skiprows=2, usecols=(5,))

data_2 = np.column_stack((rRel,y1r, y2r, y3r, y4r))
data_2 = data_2[np.argsort(data_2[:, 0])]


data = data_1
rRel = data[:,0]
y1r  = data[:,1]
y2r  = data[:,2]
y3r  = data[:,3]
y4r  = data[:,4]

ax1yr.plot(rRel, y1r, color='k')
ax2yr.plot(rRel, y2r, color='k')
ax3yr.plot(rRel, y3r, color='k')
ax4yr.plot(rRel, y4r, color='k')


data = data_2
rRel = data[:,0]
y1r  = data[:,1]
y2r  = data[:,2]
y3r  = data[:,3]
y4r  = data[:,4]

ax1yr.plot(rRel, y1r, ":", color='r')
ax2yr.plot(rRel, y2r, ":", color='r')
ax3yr.plot(rRel, y3r, ":", color='r')
ax4yr.plot(rRel, y4r, ":", color='r')


ax1yr.set_xlim(xmin, xmax)

ax1yr.set_ylabel("$y_{1, R}$")
ax2yr.set_ylabel("$y_{2, R}$")
ax3yr.set_ylabel("$y_{3, R}$")
ax4yr.set_ylabel("$y_{4, R}$")

ax1yr.set_ylim(-20, 20)
ax1yr.axhline(y=0.002,xmin=xmin,xmax=xmax, c="blue",linewidth = 0.5,zorder = 0)

ax4yr.set_xlabel("r/R")

print PSfile
plt.savefig(PSfile)
plt.close('all')
