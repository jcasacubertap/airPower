#!/usr/bin/env bash
# Build MUMPS 5.3.5 and its MATLAB MEX for DeHNSSo — macOS, Apple Silicon.
#
#   ./install_mumps_macos.sh [build_dir] [matlab_root]
#   ./install_mumps_macos.sh --ilp64 [build_dir] [matlab_root]
#
# Default (LP64): 32-bit integers, Accelerate BLAS, Homebrew METIS. This is the
# build DeHNSSo ships with and is enough for every case in the paper.
#
# --ilp64: 64-bit integers, MATLAB's own BLAS and a private 64-bit METIS, as the
# Linux build. Needed only for a Jacobian with more than 2^31 nonzeros (about
# 100 GB to assemble in MATLAB, so a workstation, not a laptop); the gateway is
# then patched for 64-bit counts as on Linux.
#
# Needs:  brew install gcc metis   and the Xcode command line tools.
# Defaults: ~/mumps_build and /Applications/MATLAB_R2025b.app
#
# Only the double-real arithmetic is built, which is all DeHNSSo uses.
set -euo pipefail

ILP64=0
if [ "${1:-}" = "--ilp64" ]; then ILP64=1; shift; fi

BUILD_DIR="${1:-$HOME/mumps_build}"
MATLAB_ROOT="${2:-/Applications/MATLAB_R2025b.app}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER=5.3.5
SRC="$BUILD_DIR/MUMPS_$VER"
METIS64DIR="$BUILD_DIR/metis64"
MEXEXT=mexmaca64
# where Homebrew's gfortran runtime lives (the version moves with every gcc bump)
GFORTLIBDIR="$(dirname "$(gfortran -print-file-name=libgfortran.dylib)")"

echo "==> build dir : $BUILD_DIR"
echo "==> MATLAB    : $MATLAB_ROOT"
echo "==> integers  : $([ $ILP64 = 1 ] && echo '64-bit (ILP64)' || echo '32-bit (LP64)')"

for t in gfortran clang make curl; do
    command -v $t >/dev/null || { echo "ERROR: $t not found."; exit 1; }
done
[ -x "$MATLAB_ROOT/bin/mex" ] || { echo "ERROR: no mex at $MATLAB_ROOT/bin/mex"; exit 1; }
if [ $ILP64 = 0 ]; then
    [ -f /opt/homebrew/opt/metis/lib/libmetis.dylib ] || \
        { echo "ERROR: METIS not found. Run: brew install metis"; exit 1; }
fi

mkdir -p "$BUILD_DIR"
if [ ! -d "$SRC" ]; then
    echo "==> downloading MUMPS $VER (public, no registration)"
    curl -fL -o "$BUILD_DIR/MUMPS_$VER.tar.gz" "https://mumps-solver.org/MUMPS_$VER.tar.gz"
    tar xf "$BUILD_DIR/MUMPS_$VER.tar.gz" -C "$BUILD_DIR"
fi

if [ $ILP64 = 1 ] && [ ! -f "$METIS64DIR/lib/libmetis.a" ]; then
    echo "==> building METIS 5.1.0 with 64-bit indices"
    # Homebrew's libmetis is IDXTYPEWIDTH 32 and crashes an ILP64 MUMPS
    cd "$BUILD_DIR"
    [ -d metis-5.1.0 ] || {
        curl -fL -o metis-5.1.0.tar.gz \
          "https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz"
        tar xf metis-5.1.0.tar.gz; }
    cd metis-5.1.0
    sed -i '' 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/' include/metis.h
    grep -q '#define IDXTYPEWIDTH 64' include/metis.h || {
        echo "ERROR: could not set METIS IDXTYPEWIDTH to 64."; exit 1; }
    # cmake directly, not METIS's `make config`: 5.1.0 predates CMake 3.5 (removed
    # in CMake 4) and its thread-local-storage probe cannot find its own source
    # file under the newer try_compile, so both are supplied here.
    cmake -S . -B build/i64 -DGKLIB_PATH="$PWD/GKlib" \
          -DCMAKE_INSTALL_PREFIX="$METIS64DIR" -DCMAKE_C_COMPILER=clang \
          -DCMAKE_C_FLAGS="-fPIC -Wno-implicit-function-declaration" \
          -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DHAVE_THREADLOCALSTORAGE=1
    cmake --build build/i64 -j4 && cmake --install build/i64
    [ -f "$METIS64DIR/lib/libmetis.a" ] || {
        echo "ERROR: METIS built but libmetis.a is missing."; exit 1; }
    cd "$SRC"
fi

echo "==> installing Makefile.inc"
# one template per platform; --ilp64 rewrites the six lines that differ
cp "$HERE/Makefile.inc.macos-arm64" "$SRC/Makefile.inc"
if [ $ILP64 = 1 ]; then
    sed -i '' \
      -e "s#^LMETISDIR = .*#LMETISDIR = $METIS64DIR/lib#" \
      -e "s#^IMETIS    = .*#IMETIS    = -I$METIS64DIR/include#" \
      -e "s#^LAPACK  = .*#LAPACK  = -L$MATLAB_ROOT/bin/maca64 -lmwlapack -lmwblas#" \
      -e "s#^LIBBLAS = .*#LIBBLAS = -L$MATLAB_ROOT/bin/maca64 -lmwlapack -lmwblas#" \
      -e "s#^CDEFS   = .*#CDEFS   = -DAdd_ -DINTSIZE64#" \
      -e "s#^OPTF    = .*#OPTF    = -O -fallow-argument-mismatch -fdefault-integer-8#" \
      "$SRC/Makefile.inc"
fi

echo "==> building the library (make d)"
( cd "$SRC" && make d )

if [ $ILP64 = 1 ]; then
    echo "==> patching mumpsmex.c for ILP64 and 64-bit nonzero counts (see mumpsmex_ilp64.patch)"
    G="$SRC/MATLAB/mumpsmex.c"
    grep -q 'MUMPS_INT);' "$G" || sed -i '' \
      -e 's#MYMALLOC((dmumps_par->jcn),(int)netrue,int);#MYMALLOC((dmumps_par->jcn),(int)netrue,MUMPS_INT);#' \
      -e 's#MYMALLOC((dmumps_par->irn),(int)netrue,int);#MYMALLOC((dmumps_par->irn),(int)netrue,MUMPS_INT);#' \
      -e 's#(dmumps_par->perm_in),int,#(dmumps_par->perm_in),MUMPS_INT,#' \
      -e 's#(dmumps_par->listvar_schur),int,#(dmumps_par->listvar_schur),MUMPS_INT,#' \
      "$G"
    grep -q 'mwSize i,j,pos;' "$G" || sed -i '' \
      -e 's#^  int i,j,pos;#  mwSize i,j,pos;#' \
      -e 's#MYMALLOC((dmumps_par->jcn),(int)netrue,MUMPS_INT);#MYMALLOC((dmumps_par->jcn),(size_t)netrue,MUMPS_INT);#' \
      -e 's#MYMALLOC((dmumps_par->irn),(int)netrue,MUMPS_INT);#MYMALLOC((dmumps_par->irn),(size_t)netrue,MUMPS_INT);#' \
      -e 's#MYMALLOC((dmumps_par->a),(int)netrue,double2);#MYMALLOC((dmumps_par->a),(size_t)netrue,double2);#' \
      -e 's#dmumps_par->nz_alloc = (int)netrue;#dmumps_par->nz_alloc = (MUMPS_INT)netrue;#' \
      -e 's#dmumps_par->nz = (int)ne;#dmumps_par->nz = (MUMPS_INT)ne;#' \
      -e 's#(dmumps_par->irn)\[j\] = ((int)irn_in\[j\])+1;#(dmumps_par->irn)[j] = ((MUMPS_INT)irn_in[j])+1;#' \
      -e 's#(dmumps_par->irn)\[pos\] = (int)irn_in\[j\]+1;#(dmumps_par->irn)[pos] = (MUMPS_INT)irn_in[j]+1;#' \
      "$G"
    for pat in 'mwSize i,j,pos;' '(size_t)netrue' '(MUMPS_INT)netrue' '(MUMPS_INT)ne;' \
               '(dmumps_par->perm_in),MUMPS_INT,' '(dmumps_par->listvar_schur),MUMPS_INT,'; do
        grep -qF "$pat" "$G" || { echo "ERROR: gateway patch did not apply ($pat missing)."; exit 1; }
    done
fi

echo "==> installing MATLAB/make.inc"
sed -e "s#\$(HOME)/mumps_build/MUMPS_$VER#$SRC#" \
    -e "s#\$(GFORTLIBDIR)#$GFORTLIBDIR#g" \
    -e "s#/Applications/MATLAB_R2025b.app#$MATLAB_ROOT#" \
    "$HERE/make.inc.matlab.macos-arm64" > "$SRC/MATLAB/make.inc"
if [ $ILP64 = 1 ]; then
    sed -i '' \
      -e "s#^LMETISDIR  = .*#LMETISDIR  = $METIS64DIR/lib#" \
      -e "s#^LIBBLAS = .*#LIBBLAS = -L$MATLAB_ROOT/bin/maca64 -lmwlapack -lmwblas#" \
      -e "s#^OPTC = .*#OPTC = -O -DINTSIZE64#" \
      "$SRC/MATLAB/make.inc"
fi

echo "==> building the MEX (d only)"
( cd "$SRC/MATLAB" && rm -f dmumpsmex.$MEXEXT dmumpsmex.stamp && make d )

MEXFILE="$SRC/MATLAB/dmumpsmex.$MEXEXT"
[ -f "$MEXFILE" ] || { echo "ERROR: mex was not produced."; exit 1; }

echo "==> checking what it linked"
otool -L "$MEXFILE" | grep -qi openblas && { echo "ERROR: linked OpenBLAS."; exit 1; }
if [ $ILP64 = 1 ]; then
    otool -L "$MEXFILE" | grep -qi libmwblas || { echo "ERROR: MATLAB BLAS not linked."; exit 1; }
    echo "    MATLAB BLAS OK, no OpenBLAS"
else
    otool -L "$MEXFILE" | grep -qi Accelerate || { echo "ERROR: Accelerate not linked."; exit 1; }
    echo "    Accelerate OK, no OpenBLAS"
fi

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
