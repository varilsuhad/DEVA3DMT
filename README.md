# DEVA3DMT
3D Magnetotelluric inversion code

**Author** : Deniz Varilsuha

**Email** : deniz.varilsuha@itu.edu.tr

This work 

## License

DEVA3DMT is released under a **non-commercial research license**.
Commercial use is prohibited without written permission from the author.
See the [LICENSE](LICENSE) file for full terms.



# DEVA3DMT (light version)

DEVA3DMT is a 3-D magnetotelluric (MT) inversion framework based on
a hybrid finite-element / finite-difference forward engine on
hexahedral meshes and GPU-accelerated iterative solvers.  
The original, full-featured code is described in:

- Varılsüha (2020), *Geophysics* 85(5): E191–E205 :contentReference[oaicite:0]{index=0}  
- Varılsüha (2024), *Computers & Geosciences* 188: 105614 :contentReference[oaicite:1]{index=1}  

This public repository contains a **lightweight, single-GPU version**
of DEVA3DMT intended for research reproducibility and educational use.
Several advanced features of the full research code have been removed
or simplified to make the package easier to compile, run, and inspect.

If you are interested in **production-level 3-D MT inversion with
all features enabled** (multi-GPU, mesh decoupling, GUI, etc.),
please contact the author (see below).

---

## Main features

### Hybrid FE–FD forward modelling

- Solves the 3-D MT forward problem using the **ungauged vector and
  scalar potential formulation** (A–φ) on structured hexahedral meshes. :contentReference[oaicite:2]{index=2}  
- **Hybrid numerical engine (HYB)**:
  - Distorted hexahedra near topography are treated with
    edge-based finite elements.
  - Undistorted rectangular blocks are treated with a finite-difference
    stencil assembled in an FE-like way. :contentReference[oaicite:3]{index=3}  
  - FE and FD contributions are combined into a single global sparse
    matrix for each polarization. :contentReference[oaicite:4]{index=4}  
- The hybrid approach is **as accurate as pure FE**, but with
  substantially fewer nonzeros and lower memory footprint, enabling
  larger models on a single GPU. :contentReference[oaicite:5]{index=5}  

### Inversion capabilities

The inversion engine follows the methodology in Varılsüha (2020). :contentReference[oaicite:6]{index=6}  

- Minimizes a composite objective function combining
  data misfit and model roughness, with optional distortion terms. :contentReference[oaicite:7]{index=7}  
- Uses **L-BFGS** as the outer optimization algorithm for log-conductivity
  and (optionally) distortion parameters, with an efficient line search. :contentReference[oaicite:8]{index=8}  
- Forward and pseudo-forward problems are solved iteratively
  (Krylov methods) and accelerated on a single GPU. :contentReference[oaicite:9]{index=9}  

#### Supported data types

The inversion can work with several MT data representations: :contentReference[oaicite:10]{index=10}  

- Impedance tensor **Z**  
- Magnetic transfer functions (tipper) **W**  
- Phase tensor **Φ** (impedance-based, distortion-free)  
- Phase vector **Ψ** (MTF-based, distortion-free)  
- Intersite tensors (**Q, T, M**) and their corresponding phase tensors  

Selected combinations of these data types can be inverted jointly.

#### Distortion handling

- Optionally estimates a **galvanic distortion matrix** acting on
  impedance and MTF data. :contentReference[oaicite:11]{index=11}  
- Distortion parameters can be configured as real/complex and
  frequency-independent/frequency-dependent within the inversion
  framework. :contentReference[oaicite:12]{index=12}  

---

## What is included in this public *light* version?

This repository contains a **reduced but functional** version of
DEVA3DMT, sufficient to:

- Build structured hexahedral meshes with topography, using the
  hybrid FE–FD forward solver described in Varılsüha (2024). :contentReference[oaicite:13]{index=13}  
- Run 3-D MT forward modelling for test models (e.g., double-brick,
  hill, two-mountain).
- Perform 3-D inversion on a **single mesh shared by all frequencies**,
  using the L-BFGS-based inversion engine.
- Use a **single NVIDIA GPU** for forward and pseudo-forward solutions.

The emphasis is on clarity and reproducibility of the core algorithms,
rather than maximum performance.

### Deliberately removed / simplified features

Compared to the full research code used in the papers, this public
version does **not** include:

- **Mesh decoupling / moving footprint inversion**  
  - In the full code, different frequency groups are solved on
    different forward meshes and mapped onto an inversion mesh, giving
    large speedups. :contentReference[oaicite:14]{index=14}  
  - Here, a *single* forward/inversion mesh is used for all frequencies.
- **Multi-GPU support**  
  - The research version distributes frequency groups across multiple
    GPUs; this version targets a **single GPU** for simplicity. :contentReference[oaicite:15]{index=15}  
- **Graphical user interface (GUI)**  
  - The private codebase contains a GUI to build models and inversion
    setups. The GUI is **not** included here; input is configured via
    MATLAB scripts and simple ASCII files.
- Some internal utilities, experimental features, and convenience
  scripts that are not essential for reproducing the published results.

If you need these capabilities (e.g., large multi-GPU inversions for
real survey data), please contact the author to discuss collaboration
or access to the full research version.

---

## Repository layout

*(Folder names may differ slightly; adjust as needed for your tree.)*

- `src/forward/` – Hybrid FE–FD forward modelling routines (A–φ
  formulation, assembly, GPU iterative solver).
- `src/inversion/` – Inversion driver, data misfit/regularization,
  L-BFGS update, sensitivity computations.
- `src/util/` – Mesh handling, I/O helpers, plotting utilities, etc.
- `examples/` – Example models and scripts:
  - `examples/double_brick/`
  - `examples/hill/`
  - `examples/two_mountain/`
- `LICENSE` – Non-commercial research license.
- `README.md` – This file.

---

## Requirements

- **MATLAB** (R20xx or later recommended) with:
  - Parallel Computing Toolbox (for GPU arrays / `gpuArray`).
- **NVIDIA GPU** with CUDA support  
  - Single device is used by this public version.
- Sufficient GPU memory to hold the global matrix and preconditioner
  for your mesh.

CPU-only runs may be possible but are not the intended or tested
configuration and will be significantly slower.

---

## Getting started

1. **Clone the repository**

   ```bash
   git clone https://github.com/<your-username>/DEVA3DMT-light.git

