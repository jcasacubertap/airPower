# Dependencies

Everything you need before running DeHNSSo, in the order you need it.

## Contents

- [1. Prerequisites](#1-prerequisites)
- [2. The `lapack_band` MEX](#2-the-lapack_band-mex)
- [3. MUMPS, for the Newton solver](#3-mumps-for-the-newton-solver)
    - [macOS (Apple Silicon)](#macos-apple-silicon)
    - [Linux (x86-64)](#linux-x86-64)
    - [What has been tested](#what-has-been-tested)

---

## 1. Prerequisites

- MATLAB R2020b or later. No toolbox is required; the Parallel Computing Toolbox is used only when `Opt.parfor = 'on' (auto-detected).
- A base-flow input file `DeHNSSo_input_<case>.mat` holding the `StabGrid` struct.
  Generate one with the scripts in `gridgen/benchmark/`.

Two optional pieces are built once per machine, and only if you need them:

| | needed for | section |
|---|---|---|
| `lapack_band` MEX | `Opt.lu_mode = 'lapack_band'`, the fastest back-end for the compressible operator | [2](#2-the-lapack_band-mex) |
| MUMPS | `Opt.solver = 'newton'`. Nothing else uses it | [3](#3-mumps-for-the-newton-solver) |

Picard runs need neither.

---

## 2. The `lapack_band` MEX

`lu_mode = 'lapack_band'` calls a MEX wrapper around LAPACK's banded LU
(`zgbsv`/`dgbsv`). The MEX binary is platform-specific and is **not** shipped
with the repo — each user builds it once.

```matlab
>> run('src/linsolve/lapack_band/mexbuild_dgbsv.m')
```

Requirements:

- MATLAB R2018a or later (uses the `-R2018a` interleaved-complex API).
- A C compiler configured for MEX (`mex -setup` once if you've never used MEX).
  - macOS: Xcode Command Line Tools.
  - Linux: `gcc`.
  - Windows: MSVC / MinGW.

If `lu_mode = 'lapack_band'` is selected but the MEX is missing, `solve_modes`
errors with a message pointing to `mexbuild_dgbsv.m`. Fall back to
`lu_mode = 'auto'` if you cannot build (e.g. no compiler available).

---

## 3. MUMPS, for the Newton solver

Newton factorizes the Jacobian with MUMPS. It is a separate dependency, and it is the only
part of DeHNSSo that needs a Fortran toolchain. macOS (Apple Silicon) and Linux (x86-64) are
both covered, with one installer each.

### macOS (Apple Silicon)

Run the installer:

```bash
src/external/mumps_5.3.5/install/install_mumps_macos.sh     # ~/mumps_build by default
```

It downloads MUMPS 5.3.5, builds the library and the MEX, checks what it linked, and solves
a test system. At the end it prints the two lines to add to your caller. If any step fails it
stops, so you never get a build it did not check. You need `brew install gcc metis` and the
Xcode command line tools.

### Linux (x86-64)

```bash
src/external/mumps_5.3.5/install/install_mumps_linux.sh     # ~/mumps_build by default
```

You need `gcc`, `gfortran`, `cmake`, `make` and `curl`. The script finds MATLAB on your
`PATH`, or you can pass its root as the second argument.

Linux is built **ILP64** (64-bit integers).

The templates are in `install/`, and their comments explain each setting.

### What has been tested

| | macOS (Apple Silicon) | Linux (x86-64) |
|---|---|---|
| integers / BLAS | LP64, Accelerate | ILP64, MATLAB BLAS |
| installer self-test (200k system) | 0.41 s, err 2.5e-08 | 0.74 s, err 2.4e-08 |
| Newton, Chebyshev grid (AMD) | `max A(1,0) = 2.641824e-02` | `max A(1,0) = 2.641824e-02` |
| Newton, FD grid (METIS) | `max A(1,0) = 2.404570e-02` | `max A(1,0) = 2.404570e-02` |

### Other architectures

Please contact the authors.

