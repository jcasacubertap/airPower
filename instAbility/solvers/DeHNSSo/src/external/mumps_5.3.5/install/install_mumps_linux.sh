#!/usr/bin/env bash
# Build MUMPS 5.3.5 and its MATLAB MEX for DeHNSSo — Linux x86-64.
#
#   ./install_mumps_linux.sh [build_dir] [matlab_root]
#
# ILP64 (64-bit integers), MATLAB's BLAS, and METIS 5 built here with 64-bit
# indices. ILP64 is not optional: the factors of a 1M-DOF Jacobian hold more
# than 2^31 entries, so a 32-bit build overflows in analysis and segfaults.
#
# Needs: gcc, gfortran, cmake, make, curl, and MATLAB.
# Defaults: ~/mumps_build, and the MATLAB found on PATH.
#
# Only the double-real arithmetic is built, which is all DeHNSSo uses.
set -euo pipefail

BUILD_DIR="${1:-$HOME/mumps_build}"
MATLAB_ROOT="${2:-$(dirname "$(dirname "$(readlink -f "$(command -v matlab)")")")}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER=5.3.5
SRC="$BUILD_DIR/MUMPS_$VER"
METIS64DIR="$BUILD_DIR/metis64"
MEXEXT=mexa64

echo "==> build dir : $BUILD_DIR"
echo "==> MATLAB    : $MATLAB_ROOT"

for t in gfortran gcc cmake make curl; do
    command -v $t >/dev/null || { echo "ERROR: $t not found."; exit 1; }
done
[ -x "$MATLAB_ROOT/bin/mex" ] || { echo "ERROR: no mex at $MATLAB_ROOT/bin/mex"; exit 1; }

mkdir -p "$BUILD_DIR"
if [ ! -d "$SRC" ]; then
    echo "==> downloading MUMPS $VER (public, no registration)"
    curl -fL -o "$BUILD_DIR/MUMPS_$VER.tar.gz" "https://mumps-solver.org/MUMPS_$VER.tar.gz"
    tar xf "$BUILD_DIR/MUMPS_$VER.tar.gz" -C "$BUILD_DIR"
fi

if [ ! -f "$METIS64DIR/lib/libmetis.a" ]; then
    echo "==> building METIS 5.1.0 with 64-bit indices"
    # the system libmetis is 32-bit and crashes an ILP64 MUMPS
    cd "$BUILD_DIR"
    [ -d metis-5.1.0 ] || {
        curl -fL -o metis-5.1.0.tar.gz \
          "https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz"
        tar xf metis-5.1.0.tar.gz; }
    cd metis-5.1.0
    # 5.1.0 sets the index width in the header, not on the command line
    sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/' include/metis.h
    grep -q '#define IDXTYPEWIDTH 64' include/metis.h || {
        echo "ERROR: could not set METIS IDXTYPEWIDTH to 64."; exit 1; }
    # do NOT pass shared=: its Makefile tests 'ifneq ($(shared), not-set)', so
    # even shared=0 builds a .so. Omitting it gives the static lib we need.
    make config prefix="$METIS64DIR" cc=gcc CFLAGS="-fPIC"
    make -j4 && make install
    [ -f "$METIS64DIR/lib/libmetis.a" ] || {
        echo "ERROR: METIS built but libmetis.a is missing (a .so is not usable here)."; exit 1; }
    cd "$SRC"
fi

echo "==> installing Makefile.inc"
sed -e "s#\$(METIS64DIR)#$METIS64DIR#g" \
    -e "s#\$(MATLABROOT)#$MATLAB_ROOT#g" \
    "$HERE/Makefile.inc.linux-ilp64" > "$SRC/Makefile.inc"

echo "==> building the library (make d)"
( cd "$SRC" && make d )

echo "==> patching mumpsmex.c for ILP64 (see mumpsmex_ilp64.patch)"
# the stock gateway allocates 4-byte int where the ILP64 library writes 8-byte
# MUMPS_INT: a heap overflow that surfaces as "free(): invalid next size"
G="$SRC/MATLAB/mumpsmex.c"
grep -q 'MUMPS_INT);' "$G" || {
    sed -i.orig \
      -e 's#MYMALLOC((dmumps_par->jcn),(int)netrue,int);#MYMALLOC((dmumps_par->jcn),(int)netrue,MUMPS_INT);#' \
      -e 's#MYMALLOC((dmumps_par->irn),(int)netrue,int);#MYMALLOC((dmumps_par->irn),(int)netrue,MUMPS_INT);#' \
      -e 's#(dmumps_par->perm_in),int,#(dmumps_par->perm_in),MUMPS_INT,#' \
      -e 's#(dmumps_par->listvar_schur),int,#(dmumps_par->listvar_schur),MUMPS_INT,#' \
      "$G"
    grep -q 'MUMPS_INT' "$G" || { echo "ERROR: gateway patch did not apply."; exit 1; }
}
# the stock gateway also counts nonzeros in C int: array lengths, nz and the copy
# loop counters are cast to 32 bits, so a Jacobian with more than 2^31 nonzeros
# fails with "Malloc failed in mumpsmex.c". Widen them to the 64-bit types the
# ILP64 library and MATLAB already use.
grep -q 'mwSize i,j,pos;' "$G" || {
    sed -i \
      -e 's#^  int i,j,pos;#  mwSize i,j,pos;#' \
      -e 's#MYMALLOC((dmumps_par->jcn),(int)netrue,MUMPS_INT);#MYMALLOC((dmumps_par->jcn),(size_t)netrue,MUMPS_INT);#' \
      -e 's#MYMALLOC((dmumps_par->irn),(int)netrue,MUMPS_INT);#MYMALLOC((dmumps_par->irn),(size_t)netrue,MUMPS_INT);#' \
      -e 's#MYMALLOC((dmumps_par->a),(int)netrue,double2);#MYMALLOC((dmumps_par->a),(size_t)netrue,double2);#' \
      -e 's#dmumps_par->nz_alloc = (int)netrue;#dmumps_par->nz_alloc = (MUMPS_INT)netrue;#' \
      -e 's#dmumps_par->nz = (int)ne;#dmumps_par->nz = (MUMPS_INT)ne;#' \
      -e 's#(dmumps_par->irn)\[j\] = ((int)irn_in\[j\])+1;#(dmumps_par->irn)[j] = ((MUMPS_INT)irn_in[j])+1;#' \
      -e 's#(dmumps_par->irn)\[pos\] = (int)irn_in\[j\]+1;#(dmumps_par->irn)[pos] = (MUMPS_INT)irn_in[j]+1;#' \
      "$G"
    grep -q 'mwSize i,j,pos;' "$G" || { echo "ERROR: 64-bit nnz patch did not apply."; exit 1; }
    grep -q '(size_t)netrue' "$G" || { echo "ERROR: 64-bit nnz patch did not apply."; exit 1; }
}

echo "==> installing MATLAB/make.inc"
sed -e "s#\$(MUMPSDIR)#$SRC#g" \
    -e "s#\$(MATLABROOT)#$MATLAB_ROOT#g" \
    -e "s#\$(METIS64DIR)#$METIS64DIR#g" \
    "$HERE/make.inc.matlab.linux-ilp64" > "$SRC/MATLAB/make.inc"

echo "==> building the MEX (d only)"
( cd "$SRC/MATLAB" && rm -f dmumpsmex.$MEXEXT dmumpsmex.stamp && make d )

MEXFILE="$SRC/MATLAB/dmumpsmex.$MEXEXT"
[ -f "$MEXFILE" ] || { echo "ERROR: mex was not produced."; exit 1; }

echo "==> checking what it linked"
ldd "$MEXFILE" | grep -qi libmwblas || { echo "ERROR: MATLAB BLAS not linked."; exit 1; }
ldd "$MEXFILE" | grep -qi libmwmetis && { echo "ERROR: bound MATLAB's METIS 4; the static link failed."; exit 1; }
echo "    MATLAB BLAS OK, METIS is static"

echo "==> solving a test system in MATLAB"
TESTM="$SRC/MATLAB/dehnsso_mumps_selftest.m"
cat > "$TESTM" <<'TESTEOF'
n = 200000;  e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);   b = A*ones(n,1);
id = initmumps;  id.SYM = 0;  id = dmumps(id);
id.JOB = 6;  id.RHS = b;
t = tic;  id = dmumps(id, A);  dt = toc(t);
err = max(abs(id.SOL - 1));
id.JOB = -2;  id = dmumps(id);          % capture the output or the instance leaks
fprintf('    solved %dx%d in %.2f s, max err %.2e\n', n, n, dt, err);
if err > 1e-6
    error('MUMPS self-test: solution is wrong (err %.2e)', err);
end
TESTEOF
OMP_NUM_THREADS=1 "$MATLAB_ROOT/bin/matlab" -batch \
    "addpath('$SRC/MATLAB'); run('$TESTM')" || { echo "ERROR: test solve failed."; exit 1; }
rm -f "$TESTM"

cat <<MSG

Done. Add this to your caller, AFTER any genpath of src/ (genpath would
otherwise shadow it with src/external/mumps_5.3.5/dmumps.m):

    addpath('$SRC/MATLAB', '-begin');
    fprintf('dmumps resolves to: %s\n', which('dmumps'));

and always start MATLAB with OMP_NUM_THREADS=1:

    OMP_NUM_THREADS=1 matlab -batch "maxNumCompThreads(1); my_caller"
MSG
