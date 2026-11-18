#!/bin/bash
# DeStreak — sigmoidal forward-facing STEP sweep (strip-like), varying width.
# Steep C2 rise (front exponent p) + plateau/gentle fall (back exponent q).
set -u
PG=/Users/jcasacubertap/OpenFOAM/cases/variants/playground/testDestreaking/airPower
DEST=/Users/jcasacubertap/OpenFOAM/cases/canonical/airPower/instAbility/TrOY/OptimalDecayNumerics/DeStreak
INP=$PG/inputs.jl
MID=$PG/PreProcessing/modules/directFlatPlateModule/postProcessing/midPlane.csv
export AIRPOWER_MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab GKSwstype=nul
JULIA=/Users/jcasacubertap/.juliaup/bin/julia
sf(){ perl -pi -e "if(/$1/){ s/(\\b$2\\s*=\\s*)[^,\\s]+/\${1}$3/ }" "$INP"; }

# id   A     xStart  xPeak   xEnd    p  q      (all shape=:sigmoidal, A>0 protrusion step)
CFG=(
  "sN  1E-3  0.150   0.158   0.190   6  3"
  "sW  1E-3  0.150   0.158   0.250   6  3"
)
echo "===== SIG step sweep start: $(date) ====="
for row in "${CFG[@]}"; do
  read -r id A xs xp xe p q <<< "$row"
  if [ -s "$DEST/bumps/bf_${id}.csv" ]; then echo "[$id] present — skip"; continue; fi
  echo ""; echo "----- [$id] step A=$A rise=[$xs,$xp] fall=[$xp,$xe] p=$p q=$q  $(date +%T) -----"
  sf "true to activate"     enabled  "true"
  sf "sigmoidal or :esn"    shape    ":sigmoidal"
  sf "positive=protrusion"  A        "$A"
  sf "front steepness"      p        "$p"
  sf "back steepness"       q        "$q"
  sf "m\\]  bump height"    A        "$A"
  sf "start of bump"        xStart   "$xs"
  sf "peak location"        xPeak    "$xp"
  sf "end of bump"          xEnd     "$xe"
  grep -nE "sigmoidal or :esn|front steepness|back steepness|start of bump|peak location|end of bump|bump height" "$INP"
  t0=$(date +%s)
  ( cd "$PG" && caffeinate -imsu "$JULIA" run.jl PreProcessing DirectFlatPlate all ) > "$DEST/dfp_${id}.log" 2>&1
  echo "[$id] exit=$? elapsed=$(( $(date +%s)-t0 ))s"
  [ -f "$MID" ] && cp "$MID" "$DEST/bumps/bf_${id}.csv" && echo "[$id] saved bf_${id}.csv" || echo "[$id] ERROR no midPlane"
done
sf "sigmoidal or :esn" shape ":esn"   # restore
echo "===== SIG step sweep done: $(date) ====="
