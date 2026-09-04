# 3d_seismic_modeling

MPI-parallel Fortran suite for 3-D frequency-domain scalar wave modeling in the
ω–kx–ky–kz domain.

The code builds a centered wavenumber grid, represents medium heterogeneity by a
spectral convolution operator, solves the frequency-domain system with Intel MKL,
and reconstructs spatial and time-domain wavefields by inverse FFTs.

- **License:** MIT
- **Language:** Fortran 90
- **Repository:** https://github.com/cssmodel/3d_seismic_modeling

## Requirements

- MPI Fortran compiler (`mpif90`; OpenMPI or MPICH with GNU Fortran or Intel Fortran)
- Intel Math Kernel Library (MKL), or a compatible BLAS/LAPACK library that provides
  `cgetrf` and `cgetrs`
- Multi-core CPU workstation or cluster (GPU is not required)

## Repository contents

| File | Role |
| --- | --- |
| `3d_make_model.f90` | Read the velocity model; set sponge boundary and source coordinates |
| `3d_full_fft.f90` | 3-D FFT of velocity, sponge, and source fields |
| `3d_full_fft.dat` | Input control file for the FFT step |
| `3d_mkl_omega_kx_ky_kz.f90` | Assemble and solve the frequency-domain convolution system |
| `3d_inverse_full_fft.f90` | 3-D inverse FFT back to space |
| `snap_syn_freq_time.f90` | Frequency-to-time seismograms and snapshots |
| `salt_256.256.064.bin` | Example 256 × 256 × 64 salt velocity model |
| `3D_Seismic_Modeling_Manual_and_License.docx` | User manual |
| `LICENSE` | MIT license |

## Workflow

Compile and run the five programs in this order:

```
3d_make_model.f90
  -> 3d_full_fft.f90
  -> 3d_mkl_omega_kx_ky_kz.f90
  -> 3d_inverse_full_fft.f90
  -> snap_syn_freq_time.f90
```

### 1. Model generation (`3d_make_model.f90`)

- **Input:** `salt_256.256.064.bin`
- **Operations:** define the spatial layout, sponge absorbing boundary, and source coordinates
- **Output:**
  - `sponge.bin` — sponge coefficients
  - `fxyz.bin` — spatial source term

### 2. Spatial-to-spectral transform (`3d_full_fft.f90`)

- **Input:** `sponge.bin`, `fxyz.bin`, and the spatial velocity model
- **Operations:** 3-D FFT and reordering onto a centered wavenumber grid
- **Output:**
  - `vp_coeff.asc` — wavenumber-domain velocity coefficients
  - `damp_coeff.asc` — wavenumber-domain damping coefficients
  - `fxyz_coeff.asc` — wavenumber-domain source coefficients

### 3. Matrix assembly and frequency-domain solve (`3d_mkl_omega_kx_ky_kz.f90`)

- **Input:** `vp_coeff.asc`, `damp_coeff.asc`, `fxyz_coeff.asc`
- **Operations:** build the impedance matrix by 3-D spectral convolution and solve
  with MKL/MPI
- Default recording length: `tmax = 10 s`
- Frequency interval: `df = 1 / tmax`
- Reference minimum velocity: 1.5 km/s (adjustable)
- Sampling: typically 2–4 grid points per wavelength

### 4. Inverse FFT (`3d_inverse_full_fft.f90`)

Returns the frequency-domain solution to the spatial domain.

### 5. Time-domain snapshots (`snap_syn_freq_time.f90`)

Converts frequency-domain results into synthetic seismograms and snapshots.

## Compile and run

```bash
mpif90 -O3 3d_make_model.f90 -o make_model
mpif90 -O3 3d_full_fft.f90 -o full_fft
mpif90 -O3 3d_mkl_omega_kx_ky_kz.f90 -lmkl_rt -o mkl_omega
mpif90 -O3 3d_inverse_full_fft.f90 -o inverse_fft
mpif90 -O3 snap_syn_freq_time.f90 -o syn_freq_time

./make_model
./full_fft
mpirun -np <number_of_processes> ./mkl_omega
./inverse_fft
./syn_freq_time
```

Change `<number_of_processes>` to the number of MPI ranks available on your machine.

## Runtime files and directories

**Generated / auxiliary inputs**

- `conect.dat` — frequency configuration (`fmax`, `nfreq`, `df`)
- `vp_coeff.asc`, `damp_coeff.asc`, `fxyz_coeff.asc`

**Outputs**

- `true_wavefield/` — frequency-domain wavefields `true.<ifreq>.<ishot>`
- `snap3d/` — complex 3-D snapshots
- `real_snap3d/` — real-valued 3-D snapshots

## Example data

`salt_256.256.064.bin` is a subsampled volume for installation tests.
It is provided so that the pipeline can be run without the full SEG/EAGE
3-D salt grid. The full benchmark model used in the paper is the publicly
available SEG/EAGE salt dataset (Aminzadeh et al., 1997).

## Citation

If you use this code, please cite the accompanying manuscript:

Shin, G. W., and Shin, C. Three-dimensional scalar wave equation modeling
in the ω–kx–ky–kz domain. *Computers & Geosciences* (submitted).

## Contact

Gee Won Shin  
Department of Industrial Engineering, Seoul National University  
Email: geewonshin@snu.ac.kr  
Phone: +82-10-7357-7558
