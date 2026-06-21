# openfoamExtensions

Custom OpenFOAM C++ extensions for airPower, compiled into `$FOAM_USER_LIBBIN`.
These are **solver-side** libraries loaded at run time by the base-flow cases
(they are not Julia preprocessing scripts).

```
openfoamExtensions/
├── Allwmake                       # builds every extension below
└── fvOptions/
    └── SFDOption/                 # Selective Frequency Damping  -> libSFDOption.so
        ├── SFDOption.{C,H}
        ├── Make/{files,options}
        └── exampleFvOptions       # ready-to-copy system/fvOptions
```

## Build (once, and after editing any source)

Inside the OpenFOAM v2312+ environment (native or the Docker container):

```sh
source $WM_PROJECT_DIR/etc/bashrc      # if not already sourced
cd PreProcessing/src/openfoamExtensions
./Allwmake
```

This installs `libSFDOption.so` into `$FOAM_USER_LIBBIN` (outside the repo), so
the build leaves no artifacts in the tracked tree except `lnInclude/` and
`Make/<arch>/` under the source dir (both git-ignorable).

## SFDOption — Selective Frequency Damping

Computes the **unstable steady base flow** of an otherwise unsteady case
(e.g. vortex shedding) by adding the forcing `-chi*(U - Ubar)` to the momentum
equation, where `Ubar` is a temporally low-pass-filtered copy of the velocity
(time constant `Delta`). As `Ubar -> U`, the flow settles onto the steady NS
solution. This is the base flow that feeds the stability analysis in
`instAbility/` (DeHNSSo / TrOY).

Encapsulated formulation: Jordi, Cotter & Sherwin (IJNMF 2014); optimal
parameter theory and the Re=100 cylinder validation: Casacuberta et al. (JCP 2018).

### Using it in a case

1. **Use a transient solver.** SFD integrates a filter ODE in real time, so the
   base-flow run must use `pimpleFoam` (NOT `simpleFoam`). That is the point of
   SFD: it reaches steady states that `simpleFoam` cannot (the unstable ones).
2. Copy `exampleFvOptions` to `<case>/system/fvOptions` and set `chi`, `Delta`.
   `fvOptions` auto-loads its own `libs ("libSFDOption.so")`.
3. `Ubar` is written with `AUTO_WRITE` (restart-safe) and the convergence
   indicator is logged every step:  `SFD: ||U - Ubar||_2 = ...`.

### Inputs (the only two physical knobs)

| entry   | meaning                    | rule of thumb                |
|---------|----------------------------|------------------------------|
| `chi`   | gain `[1/s]`               | `chi   > |omega|/2`          |
| `Delta` | filter width `[s]`         | `Delta > 1/|omega|`          |

`omega` = angular frequency of the unstable mode. For a known leading eigenvalue
`mu = mu_r + i*mu_i` (`mu_i` = growth rate, `mu_r` = frequency), the
double-root (scalar-optimal) choice is
`chi* = (|mu| + mu_i)/2`, `Delta* = 2/(|mu| - mu_i)`.

> Note: the exact asymptotic decay rate of `||U - Ubar||` is set by the
> time-integration of the host solver (alpha_num), not by SFD alone; the
> converged base flow itself is scheme-independent.
