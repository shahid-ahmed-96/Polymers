########################### Code Overview #####################################
# This script compares the Vibrational Density of States (VDOS) for different
# polymer systems by plotting single-chain (sc) and crystalline (cryst) spectra.
#
# The VDOS data should already be generated using the VACF → FFT procedure
# (see the VDOS generation script in this repository).
#
# Workflow:
#   1. Load VDOS data files for both single-chain and crystalline systems.
#   2. Remove small negative values (numerical FFT artifacts).
#   3. Normalize the VDOS using the integral of the spectrum.
#   4. Limit the frequency range to a specified maximum (freq_lim).
#   5. Smooth the spectrum using a Savitzky–Golay filter to reduce noise.
#   6. Plot and save the comparison figure.
#
# -------------------------------------------------------------------------------
# Expected input files
#
# vdos_300K_sc_<polymer>.dat
# vdos_300K_cryst_<polymer>.dat
#
# Example:
#   vdos_300K_sc_sPP.dat
#   vdos_300K_cryst_sPP.dat
#
# File format:
#   Column 1 : Frequency (THz)
#   Column 2 : VDOS
#
# -------------------------------------------------------------------------------
# Output
#
# For each polymer type, the script saves a figure:
#
#   vdos_<polymer>_<freq_lim>.png
#
# Example:
#   vdos_sPP_20.png
#
# The plot compares:
#   - Black line  : Single-chain VDOS
#   - Colored line: Crystalline VDOS
#
# Polymer color convention used in this work:
#   PP systems → firebrick
#   PS systems → teal
#
# -------------------------------------------------------------------------------
# Polymers currently analyzed:
#
#   sPP : syndiotactic polypropylene
#   iPP : isotactic polypropylene
#   sPS : syndiotactic polystyrene
#   iPS : isotactic polystyrene
#
###############################################################################


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

freq_lim = 20

# Load VDOS data
def vdos(polymername, freq_lim):
    sc_data = np.loadtxt('vdos_300K_sc_%s.dat'%polymername, skiprows=1)
    cr_data = np.loadtxt('vdos_300K_cryst_%s.dat'%polymername, skiprows=1)

    # Extract frequency and VDOS
    freq_sc = sc_data[:, 0]
    freq_cr = cr_data[:, 0]
    
    sc = sc_data[:, 1]
    cr = cr_data[:, 1]

    sc[sc < 0] = 0
    cr[cr < 0] = 0

    sc /= np.trapz(sc, freq_sc)
    cr /= np.trapz(cr, freq_cr)

    mask_sc = freq_sc <= freq_lim
    mask_cr = freq_cr <= freq_lim

    freq_sc = freq_sc[mask_sc]
    freq_cr = freq_cr[mask_cr]
    sc = sc[mask_sc]
    cr = cr[mask_cr]

    sc_smooth = savgol_filter(sc, window_length=2801, polyorder=5)
    cr_smooth = savgol_filter(cr, window_length=2801, polyorder=5)

    if polymername =='iPS':
        _ = 'teal'
    elif polymername =='sPS':
        _ = 'teal'
    else:
        _ = 'firebrick'

    plt.figure()
    plt.plot(freq_sc, sc_smooth, color='black', linewidth=2, label=' ')
    plt.plot(freq_cr, cr_smooth, color=_, linewidth=2, label=' ')
    plt.xlabel('Frequency (THz)')
    plt.ylabel('VDOS (arb. units)')
    plt.xlim(0, freq_lim)
    plt.ylim(0, )
    plt.yticks([])
    plt.savefig('vdos_%s_%d.png'%(polymername, freq_lim), bbox_inches='tight', dpi=1000, transparent=True)
    plt.close()


polymers = ['sPP', 'iPP', 'sPS', 'iPS']

for p in polymers:
    vdos(p, freq_lim)