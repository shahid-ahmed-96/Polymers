########################### Code Overview #####################################
# This script computes the Vibrational Density of States (VDOS) from velocity
# autocorrelation function (VACF) data generated from LAMMPS molecular dynamics
# simulations.
#
# The workflow is:
#   1. Read simulation parameters (temperature, seeds, timestep, units) from
#      the `input.dat` file in the current directory.
#   2. Load VACF data from each seed directory.
#   3. Average the VACF over all seeds to improve statistics.
#   4. Apply a Hann window to reduce spectral leakage before performing FFT.
#   5. Compute the Fourier transform of the VACF to obtain the VDOS.
#   6. Save the positive-frequency portion of the spectrum to `vdos.dat`.
#
# -------------------------------------------------------------------------------
# Expected directory structure
#
# Current/
# │
# ├── input.dat                  # Contains simulation parameters
# │
# ├── Temp_300/
# │   ├── seed_1/
# │   │   └── raw_vacf_1.txt
# │   ├── seed_2/
# │   │   └── raw_vacf_2.txt
# │   └── ...
#
# Each `raw_vacf_seed.txt` file is the VACF output from LAMMPS.
#
# -------------------------------------------------------------------------------
# input.dat should contain entries like:
#
# units = real
# T = [300]
# seeds = [1,2,3,4]
# timestep = 0.5
#
# -------------------------------------------------------------------------------
# Output
#
# vdos.dat
#   Column 1 : Frequency (THz)
#   Column 2 : Vibrational Density of States
#
# NOTE:
# The output file should ideally be renamed according to the system being
# analyzed (e.g., `vdos_300K_sc_sPP.dat`, `vdos_300K_cryst_iPS.dat`). See
# `compare_vacf.py` in this repository for example plotting scripts.
#
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

with open('input.dat', 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith("units ="):
            units = line.split('=')[1].strip()
        elif line.startswith("T ="):
            temps_str = line.split('=')[1].strip().strip('[]')
            Temp = list(map(int, temps_str.split(',')))
        elif line.startswith("seeds ="):
            seeds_str = line.split('=')[1].strip().strip('[]')
            seeds = list(map(int, seeds_str.split(',')))
        elif line.startswith("timestep ="):
            timestep = float(line.split('=')[1].strip())

sampling = 5
if units == 'real':
    timestep = timestep*sampling*1e-3

times = np.loadtxt(f'Temp_{Temp[0]}/seed_{seeds[0]}/raw_vacf_{seeds[0]}.txt', skiprows=2)[:,0]*timestep #Everything in ps
VACF_ = np.zeros((len(times), len(seeds)))
for i in range(len(seeds)):
    vacf_ = np.loadtxt(f'Temp_300/seed_{seeds[i]}/raw_vacf_{seeds[i]}.txt', skiprows=2)[:,-1]
    VACF_[:,i] = vacf_

VACF = np.mean(VACF_, axis=1)

vacf_windowed = VACF * np.hanning(len(VACF))

# FFT
vdos_raw = np.real(np.fft.fft(vacf_windowed))
freqs = np.fft.fftfreq(len(vacf_windowed), d=timestep)  # timestep in ps → freq in THz

vdos_raw_nofilter = np.real(np.fft.fft(VACF))
freqs_nofilter = np.fft.fftfreq(len(VACF), d=timestep)

# Positive frequencies only
positive_freqs = freqs[freqs >= 0]
vdos = vdos_raw[freqs >= 0]

positive_freqs_nofilter = freqs_nofilter[freqs_nofilter >= 0]
vdos_nofilter = vdos_raw_nofilter[freqs_nofilter >= 0]

np.savetxt('vdos.dat', np.column_stack([positive_freqs, vdos]), header='Frequency (THz)   VDOS')