#3D Seismic Modeling Package

Documentation, Manual, and Licensing Guide

The 3D Seismic Modeling Package is a high-performance 3D frequency-wavenumber

domain (ω-kx-ky-kz) seismic wave propagation modeling pipeline developed in Fortran. 

It utilizes the Intel Math Kernel Library (MKL) for fast 3D Fast Fourier Transforms (FFT) and 

handles model generation, 3D forward wavefield calculation, and synthetic seismogram/snapshot extraction.

MPI-parallelized Fortran simulation package for 3D frequency-domain wavefield modeling,

utilizing Fast Fourier Transforms (FFT) and Intel MKL BLAS/LAPACK routines for matrix

factorization and solver operations.

License: MIT 

Repository: https://github.com/cssmodel/3d_seismic_modeling

Key Features

• 3D Domain Wave Modeling: Efficient simulation of 3D acoustic wave propagation in 

frequency-wavenumber domain.

• Intel MKL Integration: Optimized FFT operations using Intel MKL libraries for high

performance scientific computing.

• Modularity: Step-by-step execution pipeline separating velocity model creation,

forward domain transformations, and time-domain snapshot syntheses.

• Salt Model Compatible: Tested and calibrated with 3D Salt velocity models (e.g.,

salt_256.256.064.bin). 

Overview and Workflow

The simulation pipeline consists of five sequential Fortran programs. They must be

compiled and executed in the exact order specified below to generate model parameters,

solve the frequency-domain wave equation, and transform wavefields back to the time

domain.

3d_make_model.f90 -> 3d_full_fft.f90 -> 3d_mkl_omega_kx_ky_kz.f90 ->

3d_inverse_full_fft.f90 -> snap_syn_freq_time.f90

Step Source File

Description / Function

Step

1

3d_make_model.f90

3d_full_fft.f90

Generates the 3D velocity grid model (e.g., binary format salt model).

Performs full 3D Fast Fourier Transform on the spatial velocity/wavefield model.

Step 3d_mkl_omega_kx_ky_kz.f90 Core solver computing wave propagation in the 3D Step Source File

2

Description / Function

frequency-wavenumber domain using Intel MKL.

3

3d_inverse_full_fft.f90

snap_syn_freq_time.f90

Transforms frequency-wavenumber domain results

back to spatial domain (Inverse 3D FFT).

Extracts synthetic seismograms and time-slice

snapshots from simulated wavefield.

Program Execution Sequence and Specifications

1. Model Generation (3d_make_model.f90)
2.  
• Input: salt_256.256.064.bin (3D salt velocity model binary).
 
• Operations: Defines the spatial layout, sponge absorbing boundary conditions, and source coordinates.

• Outputs: –

sponge.bin: Spatial sponge boundary coefficient matrix. –

fxyz.bin: Spatial source term distribution matrix.

3. Spatial-to-Spectral Transformation (3d_full_fft.f90)
4.  
• Input: sponge.bin, fxyz.bin, and the spatial velocity model.
 
• Operations: Converts the spatial velocity model, sponge model, and sources into 

spatial frequencies, reordering and transforming the arrays into spectral components.

• Outputs: –

vp_coeff.asc: Wavenumber-domain velocity coefficient file. – – 

damp_coeff.asc: Wavenumber-domain damping boundary coefficient file. 

fxyz_coeff.asc: Wavenumber-domain source coefficient file. 

5. Matrix Assembly and Frequency-Domain Solver (3d_mkl_omega_kx_ky_kz.f90)
6. 
• Input: vp_coeff.asc, damp_coeff.asc, and fxyz_coeff.asc.

• Operations: Performs convolutions to construct the global impedance matrix and solves the wave equation across distributed frequencies using MKL and MPI routines. 

7.
Core Parameter Specifications 

• Time Window (tmax): Default is set to 10 s for standard recording length, but fully modifiable. 

• Frequency Interval (df): Determined directly by the time window via df = 1/tmax. 

• Maximum Frequency and Velocity: Bounded by a reference minimum velocity of 1.5 km/s (currently hardcoded, but adjustable). 

• Discretization: Standardized to 2 grid points per wavelength, with the flexibility to scale from 2 to 4 points per wavelength depending on dispersion requirements.

• Grid Spacing Formulation: Computed as dx = (2 * nkx + 1) * dkx, dynamically evaluating the dominant spatial dimension between nkx and nkz to govern grid sizing. 

8. Inverse Fourier Transformation (3d_inverse_full_fft.f90)
9. 
• Operations: Performs the 3D inverse Fast Fourier Transform to return data to the spatial domain. 

10. Time-Domain Conversion and Snapshot Generation (snap_syn_freq_time.f90)
11. 
• Operations: Converts frequency-domain results into time-domain synthetic seismograms and snapshots. 
Input and Output Files Summary 
Input Data Files 

• salt_256.256.064.bin: Example 3D salt velocity model (256 x 256 x 64). 

• 3d_full_fft.dat: Input control file for the FFT step. 

• vp_coeff.asc: P-wave velocity coefficients and grid dimensions (nxe, nye, nze, nkx,nky, nkz, spacing dx, dy, dz).

• damp_coeff.asc: Damping boundary coefficients for absorbing boundary conditions. 

• fxyz_coeff.asc: Source term coefficients distributed across spatial coordinates.

• conect.dat: Connectivity and frequency configuration file generated during execution (stores frequency parameters like 
fmax, nfreq, df). 

Output Directories and Files 

• sponge.bin: Spatial sponge boundary coefficient matrix.

• fxyz.bin: Spatial source term distribution matrix. 

• true_wavefield/: Contains computed frequency-domain wavefield outputs formatted as true... 

• snap3d/: Directory for complex 3D wavefield snapshots. 

• real_snap3d/: Directory for real-valued 3D wavefield snapshots. 

Compilation and Execution 

Prerequisites 

• MPI Fortran Compiler (mpif90, e.g., OpenMPI or MPICH with GNU Fortran or Intel Fortran), or Intel Fortran Compiler (ifx / ifort) 

• Intel Math Kernel Library (MKL) or compatible BLAS/LAPACK library providing cgetrf and cgetrs 

MPI compile commands 

mpif90 -O3 3d_make_model.f90 -o make_model 

mpif90 -O3 3d_full_fft.f90 -o full_fft 

mpif90 -O3 3d_mkl_omega_kx_ky_kz.f90 -lmkl_rt -o mkl_omega

mpif90 -O3 3d_inverse_full_fft.f90 -o inverse_fft 

mpif90 -O3 snap_syn_freq_time.f90 -o syn_freq_time 

mpirun -np  ./mkl_omega 

Intel Fortran compile commands 

ifx -O3 3d_make_model.f90 -o 3d_make_model.exe 

./3d_make_model.exe 

ifx -O3 -qmkl 3d_full_fft.f90 -o 3d_full_fft.exe 

./3d_full_fft.exe 

ifx -O3 -qmkl 3d_mkl_omega_kx_ky_kz.f90 -o 3d_mkl_omega_kx_ky_kz.exe 

./3d_mkl_omega_kx_ky_kz.exe 

ifx -O3 -qmkl 3d_inverse_full_fft.f90 -o 3d_inverse_full_fft.exe 

./3d_inverse_full_fft.exe 

ifx -O3 snap_syn_freq_time.f90 -o snap_syn_freq_time.exe 

./snap_syn_freq_time.exe 

Software License (MIT License) 

Copyright (c) 2026 3D Seismic Modeling Authors

Permission is hereby granted, free of charge, to any person obtaining a copy of this

software and associated documentation files (the Software), to deal in the Software without 

restriction, including without limitation the rights to use, copy, modify, merge, publish, 

distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom 

the Software is furnished to do so, subject to the following conditions: 

The above copyright notice and this permission notice shall be included in all copies or 

substantial portions of the Software. 

THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 

IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,

FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL

THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 

OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,

ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 

OTHER DEALINGS IN THE SOFTWARE.

Contact 

Gee Won Shin 

Department of Industrial Engineering, Seoul National University 

Email: geewonshin@snu.ac.kr 

Telephone: +82-10-7357-7558 

