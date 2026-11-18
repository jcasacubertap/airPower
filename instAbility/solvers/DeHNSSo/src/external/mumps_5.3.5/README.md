# MUMPS for DeHNSSo

MUMPS factorises the Jacobian in the Newton solver (`src/newton/mumps_solve.m`).
Nothing else in DeHNSSo uses it.

- **To install it:** `docs/dependencies.md`, or run `install/install_mumps_macos.sh`
  or `install/install_mumps_linux.sh`.
- **To run with it:** `docs/DeHNSSo_guide.md`, "Running it: MUMPS, memory and where".

This folder holds only the upstream MATLAB wrappers DeHNSSo calls (`dmumps.m`,
`initmumps.m`) and the build files in `install/`. The MEX is not in git, so your
own build must go on the path ahead of this folder.
