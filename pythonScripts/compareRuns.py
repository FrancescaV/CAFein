import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
from math import sqrt
from math import sin
from math import pow
from math import fabs
from math import log


plt.rcParams.update({'font.size': 12})

# Constants
G_cgs      =  6.67259e-8
Rsun_toCm  = 6.95660e10
Msun_toG   = 1.98855e+33
MjPerMsun    = 9.54265748e-4
secPerMin  = 60.0
minPerHours = 60.0
hoursPerDay = 24.0

secPerDay = secPerMin * minPerHours * hoursPerDay

minsPerYears = 60.0*24.0*365.242199
secPerYears = secPerMin*minsPerYears

basicAddress = "/Users/francesca/Northwestern/Research/realModels_eigenvalues_Eigenfunc_gsl_v23_zioGio_readingAllRij_withP/output/test_08262015_J0651_accR1eM14_accV1eM14_rkf45/"


basicAddress2 = "/Users/francesca/Northwestern/Research/realModels_eigenvalues_Eigenfunc_gsl_v23_zioGio_readingAllRij/output/test_09042013_J0651_accR1eM14_accV1eM14_rkf45_zioGioHope_hope/"


addressFiles_EF = basicAddress+"timescale.dat"
addressFiles_global = basicAddress+"global.dat"
plotsAddress = basicAddress+"plots/"


y1r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(2,))
y2r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(3,))
y3r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(4,))
y4r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(5,))
y5r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(6,))
y6r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(7,))
y7r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(8,))
y8r  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(9,))

y1i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(10,))
y2i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(11,))
y3i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(12,))
y4i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(13,))
y5i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(14,))
y6i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(15,))
y7i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(16,))
y8i  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(17,))

omega_R  = np.loadtxt(addressFiles_global, skiprows=1, usecols=(10,))

addressFiles_EF2 = basicAddress2+"timescale.dat"
addressFiles_global2 = basicAddress2+"global.dat"



#
#
#y1r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(2,))
#y2r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(3,))
#y3r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(4,))
#y4r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(5,))
#y5r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(6,))
#y6r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(7,))
#y7r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(8,))
#y8r_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(9,))
#
#y1i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(10,))
#y2i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(11,))
#y3i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(12,))
#y4i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(13,))
#y5i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(14,))
#y6i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(15,))
#y7i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(16,))
#y8i_2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(17,))
#
#omega_R_2  = np.loadtxt(addressFiles_global2, skiprows=1, usecols=(10,))
#
#print len(y1r_2), len(omega_R_2)

"""
    Tidal eigenfunctions at the surface
    """
PSfile = plotsAddress+"yi_real_J0651.jpeg"

plt.close('all')
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2)

ax1.plot(omega_R, y1r, label = "07/20/2015")
ax2.plot(omega_R, y2r)
ax3.plot(omega_R, y3r)
ax4.plot(omega_R, y4r)
ax5.plot(omega_R, y5r)
ax6.plot(omega_R, y6r)
ax7.plot(omega_R, y7r)
ax8.plot(omega_R, y8r)

#ax1.plot(omega_R_2, y1r_2, color="red", label = "08/04/2013")
#ax2.plot(omega_R_2, y2r_2, color="red")
#ax3.plot(omega_R_2, y3r_2, color="red")
#ax4.plot(omega_R_2, y4r_2, color="red")
#ax5.plot(omega_R_2, y5r_2, color="red")
#ax6.plot(omega_R_2, y6r_2, color="red")
#ax7.plot(omega_R_2, y7r_2, color="red")
#ax8.plot(omega_R_2, y8r_2, color="red")

plt.tight_layout()
ax1.set_ylabel(r"$y1 = \xi_r/r$")
ax2.set_ylabel(r"$y2 = \xi_h/r$")

ax3.set_ylabel(r"$y3 = \Phi'/(gr)$")
ax4.set_ylabel(r"$y4 = (1/g)d\Phi'/dr$")
ax5.set_ylabel(r"$y5 = \delta S/c_P$")
ax6.set_ylabel(r"$y6 = \delta L_R/L_R$")
ax7.set_ylabel(r"$y_7$")
ax8.set_ylabel(r"$y_8$")

ax1.set_title("Real")
ax2.set_title("Real")
ax7.set_xlabel(r"$\omega_T$")
ax8.set_xlabel(r"$\omega_T$")


plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, wspace=0.4, hspace=0.35)
ax1.set_xlim(-0.5, 0.5)

fig.savefig(PSfile)
plt.close('all')

PSfile = plotsAddress+"yi_imag_J0651.jpeg"

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(nrows=4, ncols=2)

ax1.plot(omega_R, y1i)
ax2.plot(omega_R, y2i)
ax3.plot(omega_R, y3i)
ax4.plot(omega_R, y4i)
ax5.plot(omega_R, y5i)
ax6.plot(omega_R, y6i)
ax7.plot(omega_R, y7i, label = "07/20/2015")
ax8.plot(omega_R, y8i)

#ax1.plot(omega_R_2, y1i_2, color = "red")
#ax2.plot(omega_R_2, y2i_2, color = "red")
#ax3.plot(omega_R_2, y3i_2, color = "red")
#ax4.plot(omega_R_2, y4i_2, color = "red")
#ax5.plot(omega_R_2, y5i_2, color = "red")
#ax6.plot(omega_R_2, y6i_2, color = "red")
#ax7.plot(omega_R_2, y7i_2, color = "red", label = "08/04/2013")
#ax8.plot(omega_R_2, y8i_2, color = "red")

plt.tight_layout()
ax1.set_ylabel(r"$y1 = \xi_r/r$")
ax2.set_ylabel(r"$y2 = \xi_h/r$")
ax3.set_ylabel(r"$y3 = \Phi'/(gr)$")
ax4.set_ylabel(r"$y4 = (1/g)d\Phi'/dr$")
ax5.set_ylabel(r"$y5 = \delta S/c_P$")
ax6.set_ylabel(r"$y6 = \delta L_R/L_R$")
ax7.set_ylabel(r"$y_7$")
ax8.set_ylabel(r"$y_8$")

ax1.set_title("Imaginary")
ax2.set_title("Imaginary")
ax7.set_xlabel(r"$\omega_T$")
ax8.set_xlabel(r"$\omega_T$")
ax7.legend(loc=1)

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, wspace=0.4, hspace=0.35)


fig.savefig(PSfile)
plt.close('all')


PSfile = plotsAddress+"cazzo.jpeg"
plt.close('all')
fig, ((ax1, ax2),(ax3, ax4),(ax5, ax6)) = plt.subplots(nrows=3, ncols=2)

ax1.plot(omega_R, y3r, color="k", label = "07/20/2015 model")
#ax1.plot(omega_R_2, y3r_2, "--", color = "red", label = "old model")

ax2.plot(omega_R, y4r, color="k")
#ax2.plot(omega_R_2, y4r_2, "--", color = "red")

ax3.plot(omega_R, y3i, color="k")
#ax3.plot(omega_R_2, y3i_2, "--", color = "red")

ax4.plot(omega_R, y4i, color="k")
#ax4.plot(omega_R_2, y4i_2, "--", color = "red")

#ax5.plot(omega_R_2, y3i_2, color = "red")
#ax6.plot(omega_R_2, y4i_2, color = "red")

ax1.set_ylabel(r"$y3 = \Phi'/(gr)$ real")
ax2.set_ylabel(r"$y4 = (1/g)d\Phi'/dr$ real")

ax3.set_ylabel(r"$y3 = \Phi'/(gr)$ Imaginary")
ax4.set_ylabel(r"$y4 = (1/g)d\Phi'/dr$ Imaginary")

ax5.set_ylabel(r"$y3 = \Phi'/(gr)$ Imaginary old")
ax6.set_ylabel(r"$y4 = (1/g)d\Phi'/dr$ Imaginary old")

ax3.set_xlabel(r"$\omega_T$")
ax4.set_xlabel(r"$\omega_T$")

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.9, wspace=0.5, hspace=0.5)

fig.savefig(PSfile)
plt.close('all')


"""
    timescales
    """
plt.close('all')

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

aTS  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(29,))
wTS  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(31,))
GRTS  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(32,))


#data for J0651 paper
PSfile = plotsAddress+"timescales_compareOldModelNewModel_J0651paper.jpg"
addressFiles_EF_old = "/Users/francesca/Northwestern/Research/realModels_eigenvalues_Eigenfunc_gsl_v17/output/test_0910_0912_0914_1003_1004_from2012_fullRangeLowRes/0p23Msun_HeWD/log1673/tides/totalOutput_eigenf_test_0910_0912_0914_0916_1003_1004_from2012.dat"
addressFiles_global_old = "/Users/francesca/Northwestern/Research/realModels_eigenvalues_Eigenfunc_gsl_v17/output/test_0910_0912_0914_1003_1004_from2012_fullRangeLowRes/0p23Msun_HeWD/log1673/tides/totalOutput_global_test_0910_0912_0914_0916_1003_1004_from2012.dat"


PSfile = plotsAddress+"timescales_compareOldModelNewModel.jpg"
addressFiles_EF_old = addressFiles_EF2
addressFiles_global_old = addressFiles_global2


omega_R_old  = np.loadtxt(addressFiles_global_old, skiprows=1, usecols=(10,))

aTS_old  = np.loadtxt(addressFiles_EF_old, skiprows=2, usecols=(29,))
wTS_old  = np.loadtxt(addressFiles_EF_old, skiprows=2, usecols=(31,))
GRTS_old  = np.loadtxt(addressFiles_EF_old, skiprows=2, usecols=(31,))

aTS_plus = []
aTS_minus = []
wTS_plus = []
wTS_minus = []
wR_aTS_plus = []
wR_aTS_minus = []
wR_wTS_plus = []
wR_wTS_minus = []


for i in range(0, aTS.shape[0]):
    if aTS[i] > 0:
        aTS_plus.append(log(aTS[i], 10.0))
        wR_aTS_plus.append(omega_R[i])
    if aTS[i] < 0:
        aTS_minus.append(log(fabs(aTS[i]), 10.0))
        wR_aTS_minus.append(omega_R[i])
    if wTS[i] > 0:
        wTS_plus.append(log(wTS[i], 10.0))
        wR_wTS_plus.append(omega_R[i])
    if wTS[i] < 0:
        wTS_minus.append(log(fabs(wTS[i]), 10.0))
        wR_wTS_minus.append(omega_R[i])



aTS_plus_old = []
aTS_minus_old = []
wTS_plus_old = []
wTS_minus_old = []
wR_aTS_plus_old = []
wR_aTS_minus_old = []
wR_wTS_plus_old = []
wR_wTS_minus_old = []



for i in range(0, aTS_old.shape[0]):
    if aTS_old[i] > 0:
        aTS_plus_old.append(log(aTS_old[i],10.0))
        wR_aTS_plus_old.append(omega_R_old[i])
    if aTS_old[i] < 0:
        aTS_minus_old.append(log(fabs(aTS_old[i]), 10.0))
        wR_aTS_minus_old.append(omega_R_old[i])
    if wTS_old[i] > 0:
        wTS_plus_old.append(log(wTS_old[i], 10.0))
        wR_wTS_plus_old.append(omega_R_old[i])
    if wTS_old[i] < 0:
        wTS_minus_old.append(log(fabs(wTS_old[i]), 10.0))
        wR_wTS_minus_old.append(omega_R_old[i])

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.3, hspace=0.3)

xmin = min(omega_R)
xmax = max(omega_R)
ymin = -10
ymax = 10

sizePoints = 0.5

ax1.scatter(wR_aTS_plus, aTS_plus, color='k', s = sizePoints)
ax1.scatter(wR_aTS_minus, aTS_minus, color='r', s = sizePoints)

ax2.scatter(wR_aTS_plus_old, aTS_plus_old, color='k', s = sizePoints)
ax2.scatter(wR_aTS_minus_old, aTS_minus_old, color='r', s = sizePoints)


ax3.scatter(wR_wTS_plus, wTS_plus, color='k', s = sizePoints)
ax3.scatter(wR_wTS_minus, wTS_minus, color='r', s = sizePoints)

ax4.scatter(wR_wTS_plus_old, wTS_plus_old, color='k', s = sizePoints)
ax4.scatter(wR_wTS_minus_old, wTS_minus_old, color='r', s = sizePoints)

ax1.set_xlim(xmin, xmax)
ax2.set_xlim(xmin, xmax)
ax3.set_xlim(xmin, xmax)
ax4.set_xlim(xmin, xmax)

ax1.set_ylim(ymin, ymax)
ax2.set_ylim(ymin, ymax)
ax3.set_ylim(ymin, ymax)
ax4.set_ylim(ymin, ymax)

ax1.set_title("NEW")
ax2.set_title("OLD")

ax1.set_ylabel("log($|t_{a, tid}/yr|$)")
ax3.set_ylabel(r"log($|t_{\Omega, tid}/yr|$)")

ax3.set_xlabel(r"$\omega_T$")
ax4.set_xlabel(r"$\omega_T$")

ax1.hlines(log(fabs(GRTS[0]), 10.0), xmin, xmax, linestyles=":")
ax2.hlines(log(fabs(GRTS_old[0]), 10.0), xmin, xmax, linestyles=":")
ax3.hlines(log(fabs(GRTS[0]), 10.0), xmin, xmax, linestyles=":")
ax4.hlines(log(fabs(GRTS_old[0]), 10.0), xmin, xmax, linestyles=":")

fig.savefig(PSfile)
plt.close('all')



##Not positive and negative
#"""
#    timescales
#    """
#plt.close('all')
#matplotlib.rcParams.update({'font.size': 9})
#
#
#
#
#
#PSfile = plotsAddress+"timescales_J0651.jpeg"
#
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
#
#aTS  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(29,))
#wTS  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(31,))
#GRTS  = np.loadtxt(addressFiles_EF, skiprows=2, usecols=(32,))
#
#fig, ((ax1), (ax2)) = plt.subplots(nrows=2, ncols=1)
#
#
#ax1.plot(omega_R, aTS, label = "07/20/2015")
#ax2.plot(omega_R, wTS)
#
#
#aTS2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(29,))
#wTS2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(31,))
#GRTS2  = np.loadtxt(addressFiles_EF2, skiprows=2, usecols=(32,))
#omega_R2  = np.loadtxt(addressFiles_global2, skiprows=1, usecols=(10,))
#
#
#ax1.plot(omega_R2, aTS2, "--", color="r", label = "08/04/2013")
#ax2.plot(omega_R2, wTS2, "--", color="r")
#
#ax2.legend(loc=1)
##ax1.set_xlim(0.1,0.13)
##ax2.set_xlim(0.1,0.13)
#ax1.set_ylim(-1e5, 1e5)
#ax2.set_ylim(-1e5, 1e5)
#
#
#plt.savefig(PSfile)
#plt.close()
#
