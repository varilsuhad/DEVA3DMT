
# DEVA3DMT

DEVA3DMT is a 3-D magnetotelluric (MT) inversion framework based on
a hybrid finite-element / finite-difference forward engine on
hexahedral meshes and GPU-accelerated iterative solvers.  
The original, full-featured code is described in:

- Varılsüha (2020), *Geophysics* 85(5): E191–E205 
- Varılsüha (2024), *Computers & Geosciences* 188: 105614

This public repository contains a **lightweight, single-GPU version**
of DEVA3DMT intended for research reproducibility and educational use.
Several advanced features of the full research code have been removed
or simplified to make the package easier to compile, run, and inspect.

If you are interested in **production-level 3-D MT inversion with
all features enabled** (multi-GPU, mesh decoupling, GUI, etc.),
please contact the author (see below).

**Author** : Deniz Varılsüha

**Email** : deniz.varilsuha@itu.edu.tr

---

## Main features

### Hybrid FE–FD forward modelling

- Solves the 3-D MT forward problem using the **ungauged vector and
  scalar potential formulation** (A–φ) on structured hexahedral meshes. 
- **Hybrid numerical engine (HYB)**:
  - Distorted hexahedra near topography are treated with
    edge-based finite elements.
  - Undistorted rectangular blocks are treated with a finite-difference
    stencil assembled in an FE-like way. 
  - FE and FD contributions are combined into a single global sparse
    matrix for each polarization. 
- The hybrid approach is **as accurate as pure FE**, but with
  substantially fewer nonzeros and lower memory footprint, enabling
  larger models on a single GPU.

### Inversion capabilities

The inversion engine follows the methodology in Varılsüha (2020) but in this version is updated to inexact Gauss-Newton (GN) algorithm. 

- Minimizes a composite objective function combining
  data misfit and model roughness, with optional distortion terms. 
- Uses **inexact-GN** as the outer optimization algorithm for log-conductivity
  and (optionally) distortion parameters,.
- Forward and pseudo-forward problems are solved iteratively
  (Krylov methods) and accelerated on a single GPU.

#### Supported data types

The inversion can work with several MT data representations:

- Impedance tensor **Z**  
- Magnetic transfer functions (tipper) **W**  
- Phase tensor **Φ** (impedance-based, distortion-free)  
- Phase vector **Ψ** (MTF-based, distortion-free)  
- Intersite tensors (**Q, T, M**) and their corresponding phase tensors  
- Amplitude tensor from **Z** and vector from **W**
Selected combinations of these data types can be inverted jointly.

#### Distortion handling

- Optionally estimates a **galvanic distortion matrix** acting on
  impedance and MTF data. 
- Distortion parameters can be configured as real/complex and
  frequency-independent/frequency-dependent within the inversion
  framework.

---

## What is included in this public *light* version?

This repository contains a **reduced but functional** version of
DEVA3DMT, sufficient to:

- Build structured hexahedral meshes with topography, using the
  hybrid FE–FD forward solver described in Varılsüha (2024).
- Run 3-D MT forward modelling for test models (e.g., double-brick,
  hill, two-mountain).
- Perform 3-D inversion on a **single mesh shared by all frequencies**,
  using the GN-based inversion engine.
- Use a **single NVIDIA GPU** for forward and pseudo-forward solutions.

The emphasis is on clarity and reproducibility of the core algorithms,
rather than maximum performance.

### Deliberately removed / simplified features

Compared to the full research code used in the papers, this public
version does **not** include:

- **Mesh decoupling / moving footprint inversion**  
  - In the full code, different frequency groups are solved on
    different forward meshes and mapped onto an inversion mesh, giving
    large speedups. 
  - Here, a *single* forward/inversion mesh is used for all frequencies.
- **Multi-GPU support**  
  - The research version distributes frequency groups across multiple
    GPUs; this version targets a **single GPU** for simplicity. 
- **Graphical user interface (GUI)**  
  - The private codebase contains a GUI to build models and inversion
    setups. The GUI is **not** included here; input is configured via
    MATLAB scripts and simple ASCII files.
- Some internal utilities, experimental features (*Triple Grid Design* or *Full Newton Inversion*), and convenience
  scripts that are not essential for reproducing the published results.

If you need these capabilities (e.g., large multi-GPU inversions for
real survey data), please contact the author to discuss collaboration
or access to the full research version.

---

## Repository layout


- `DEVA3DMT_start.m` – The starting point of the Matlab code.
- `DEVA3DMT_main.m` – The main function takes in the the data and performs the inversion. 
- `digermex` and `solvercuda`– mex functions called directly from Matlab for GPU operations.
- `solvercuda` - Other Matlab subfunctions that is required for the code.
- `LICENSE` – Non-commercial research license.
- `README.md` – This file.

---

## Requirements

- **MATLAB** (R2024b or later recommended) with:
  - Parallel Computing Toolbox.
  - Statistics and Machine Learning Toolbox.
- **NVIDIA GPU** with CUDA support  
  - Single device is used by this public version.
- Sufficient GPU memory to hold the global matrix and preconditioner
  for your mesh. A minimum of 24GB Vram and preferably more.

CPU-only runs may be possible but are not the intended or tested
configuration and will be significantly slower.

---

## Getting started

1. **Clone the repository**

   ```bash
   git clone https://github.com/<your-username>/DEVA3DMT-light.git

