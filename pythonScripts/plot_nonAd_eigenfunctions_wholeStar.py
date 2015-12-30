import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
from math import sqrt
from math import sin
from math import pow
from math import fabs
matplotlib.rcParams.update({'font.size': 12})

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
xmin = 0.0
xmax = 1.1

basicAddress1 = "../output/"
addressFiles1 = basicAddress1+"eigenfunctions.dat"
plotsAddress = basicAddress1+"plots/"

basicAddress2 = "../input/test_suite/non_adiabatic_eigenfunctions/Output"
addressFiles2 = basicAddress2+"/eigenfunctions_lite.dat"
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
fig, ((ax1yr, ax1yi), (ax2yr, ax2yi), (ax3yr, ax3yi), (ax4yr, ax4yi), (ax5yr, ax5yi), (ax6yr, ax6yi)) = plt.subplots(nrows=6, ncols=2, sharex = "all")
fig.set_size_inches(fig_size)
plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, wspace=0.35, hspace=0.25)


'''
........csiRel.............csiAbs (Rsun)...................y1real..................y2real..................y3real..................y4real..................y5real..................y6real..................y7real..................y8real..................y1imag..................y2imag..................y3imag..................y4imag..................y5imag..................y6imag..................y7imag..................y8imag..................Pspin(min)................rho......................g..........................P......................Gamma1................N2.....................delAd......................T (g/mol)................cP(mol/g)....................entropy................Lr.......................a/dadt_T (yr)..........e/dedt_T (yr)..........omega/dOmegadt_T(yr)..........a/dadt_GR (yr)..........argFlmk(deg)..........modFlmk..........linearityViolation (1/0)...........(1/J)dJ/dt.............JdotSpin/Jspin.............JdotOrb/Jorb.................V..........
'''
addressFiles = addressFiles1
rRel = np.loadtxt(addressFiles, skiprows=2, usecols=(0,))
y1r  = np.loadtxt(addressFiles, skiprows=2, usecols=(2,))
y2r  = np.loadtxt(addressFiles, skiprows=2, usecols=(3,))
y3r  = np.loadtxt(addressFiles, skiprows=2, usecols=(4,))
y4r  = np.loadtxt(addressFiles, skiprows=2, usecols=(5,))
y5r  = np.loadtxt(addressFiles, skiprows=2, usecols=(6,))
y6r  = np.loadtxt(addressFiles, skiprows=2, usecols=(7,))
y1i  = np.loadtxt(addressFiles, skiprows=2, usecols=(8,))
y2i  = np.loadtxt(addressFiles, skiprows=2, usecols=(9,))
y3i  = np.loadtxt(addressFiles, skiprows=2, usecols=(10,))
y4i  = np.loadtxt(addressFiles, skiprows=2, usecols=(11,))
y5i  = np.loadtxt(addressFiles, skiprows=2, usecols=(12,))
y6i  = np.loadtxt(addressFiles, skiprows=2, usecols=(13,))

data_1 = np.column_stack((rRel,y1r, y2r, y3r, y4r, y5r, y6r, y1i, y2i, y3i, y4i, y5i, y6i))
data_1 = data_1[np.argsort(data_1[:, 0])]

addressFiles = addressFiles2
rRel = np.loadtxt(addressFiles, skiprows=2, usecols=(0,))
y1r  = np.loadtxt(addressFiles, skiprows=2, usecols=(2,))
y2r  = np.loadtxt(addressFiles, skiprows=2, usecols=(3,))
y3r  = np.loadtxt(addressFiles, skiprows=2, usecols=(4,))
y4r  = np.loadtxt(addressFiles, skiprows=2, usecols=(5,))
y5r  = np.loadtxt(addressFiles, skiprows=2, usecols=(6,))
y6r  = np.loadtxt(addressFiles, skiprows=2, usecols=(7,))
y1i  = np.loadtxt(addressFiles, skiprows=2, usecols=(8,))
y2i  = np.loadtxt(addressFiles, skiprows=2, usecols=(9,))
y3i  = np.loadtxt(addressFiles, skiprows=2, usecols=(10,))
y4i  = np.loadtxt(addressFiles, skiprows=2, usecols=(11,))
y5i  = np.loadtxt(addressFiles, skiprows=2, usecols=(12,))
y6i  = np.loadtxt(addressFiles, skiprows=2, usecols=(13,))

data_2 = np.column_stack((rRel,y1r, y2r, y3r, y4r, y5r, y6r, y1i, y2i, y3i, y4i, y5i, y6i))
data_2 = data_2[np.argsort(data_2[:, 0])]


data = data_1
rRel = data[:,0]
y1r  = data[:,1]
y2r  = data[:,2]
y3r  = data[:,3]
y4r  = data[:,4]
y5r  = data[:,5]
y6r  = data[:,6]
y1i  = data[:,7]
y2i  = data[:,8]
y3i  = data[:,9]
y4i  = data[:,10]
y5i  = data[:,11]
y6i  = data[:,12]

ax1yr.plot(rRel, y1r, color='k')
ax2yr.plot(rRel, y2r, color='k')
ax3yr.plot(rRel, y3r, color='k')
ax4yr.plot(rRel, y4r, color='k')
ax5yr.plot(rRel, y5r, color='k')
ax6yr.plot(rRel, y6r, color='k')

ax1yi.plot(rRel, y1i, color='k')
ax2yi.plot(rRel, y2i, color='k')
ax3yi.plot(rRel, y3i, color='k')
ax4yi.plot(rRel, y4i, color='k')
ax5yi.plot(rRel, y5i, color='k')
ax6yi.plot(rRel, y6i, color='k')


data = data_2
rRel = data[:,0]
y1r  = data[:,1]
y2r  = data[:,2]
y3r  = data[:,3]
y4r  = data[:,4]
y5r  = data[:,5]
y6r  = data[:,6]
y1i  = data[:,7]
y2i  = data[:,8]
y3i  = data[:,9]
y4i  = data[:,10]
y5i  = data[:,11]
y6i  = data[:,12]


ax1yr.plot(rRel, y1r, ":", color='r')
ax2yr.plot(rRel, y2r, ":", color='r')
ax3yr.plot(rRel, y3r, ":", color='r')
ax4yr.plot(rRel, y4r, ":", color='r')
ax5yr.plot(rRel, y5r, ":", color='r')
ax6yr.plot(rRel, y6r, ":", color='r')

ax1yi.plot(rRel, y1i, ":", color='r')
ax2yi.plot(rRel, y2i, ":", color='r')
ax3yi.plot(rRel, y3i, ":", color='r')
ax4yi.plot(rRel, y4i, ":", color='r')
ax5yi.plot(rRel, y5i, ":", color='r')
ax6yi.plot(rRel, y6i, ":", color='r')


ax1yr.set_xlim(xmin, xmax)
ax1yi.set_xlim(xmin, xmax)

ax1yr.set_ylabel("$y_{1, R}$")
ax2yr.set_ylabel("$y_{2, R}$")
ax3yr.set_ylabel("$y_{3, R}$")
ax4yr.set_ylabel("$y_{4, R}$")
ax5yr.set_ylabel("$y_{5, R}$")
ax6yr.set_ylabel("$y_{6, R}$")

ax1yi.set_ylabel("$y_{1, I}$")
ax2yi.set_ylabel("$y_{2, I}$")
ax3yi.set_ylabel("$y_{3, I}$")
ax4yi.set_ylabel("$y_{4, I}$")
ax5yi.set_ylabel("$y_{5, I}$")
ax6yi.set_ylabel("$y_{6, I}$")

ax6yr.set_xlabel("r/R")
ax6yi.set_xlabel("r/R")

ax1yr.axhline(y=0.002,xmin=xmin,xmax=xmax, c="blue",linewidth = 0.5,zorder = 0)

print PSfile
plt.savefig(PSfile)
plt.close('all')
