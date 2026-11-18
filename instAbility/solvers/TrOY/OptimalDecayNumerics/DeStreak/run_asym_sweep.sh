#!/bin/bash
# DeStreak — asymmetric (forward-facing-step) shape sweep via DirectFlatPlate.
# Adds an epsilon (ESN skewness) field to the per-bump edits. ε<0 = forward-facing
# (steep upstream, gentle downstream). Runs in the PLAYGROUND (private _pg dirs).
set -u
PG=/Users/jcasacubertap/OpenFOAM/cases/variants/playground/testDestreaking/airPower
DEST=/Users/jcasacubertap/OpenFOAM/cases/canonical/airPower/instAbility/TrOY/OptimalDecayNumerics/DeStreak
INP=$PG/inputs.jl
MID=$PG/PreProcessing/modules/directFlatPlateModule/postProcessing/midPlane.csv
export AIRPOWER_MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab GKSwstype=nul
JULIA=/Users/jcasacubertap/.juliaup/bin/julia
set_field(){ perl -pi -e "if(/$1/){ s/(\\b$2\\s*=\\s*)[^,\\s]+/\${1}$3/ }" "$INP"; }

# id   A     R   xCenter  eps    enabled
CFG=(
  "e05  1E-3  50  0.16   -0.5   true"
  "e07  1E-3  50  0.16   -0.7   true"
)
echo "===== ASYM sweep start: $(date) ====="
for row in "${CFG[@]}"; do
  read -r id A R xc eps en <<< "$row"
  if [ -s "$DEST/bumps/bf_${id}.csv" ]; then echo "[$id] present — skip"; continue; fi
  echo ""; echo "----- [$id] A=$A R=$R xc=$xc eps=$eps  $(date +%T) -----"
  set_field "true to activate"     enabled  "$en"
  set_field "positive=protrusion"  A        "$A"
  set_field "base-to-height ratio" R        "$R"
  set_field "skewness"             epsilon  "$eps"
  set_field "m\\]  bump height"    A        "$A"
  set_field "ESN bump center"      xCenter  "$xc"
  grep -nE "skewness|bump height|ESN bump center|base-to-height" "$INP"
  t0=$(date +%s)
  ( cd "$PG" && caffeinate -imsu "$JULIA" run.jl PreProcessing DirectFlatPlate all ) > "$DEST/dfp_${id}.log" 2>&1
  echo "[$id] exit=$? elapsed=$(( $(date +%s)-t0 ))s"
  [ -f "$MID" ] && cp "$MID" "$DEST/bumps/bf_${id}.csv" && echo "[$id] saved bf_${id}.csv" || echo "[$id] ERROR no midPlane"
done
# restore epsilon to 0 (leave playground tidy)
set_field "skewness" epsilon "0.0"
echo "===== ASYM sweep done: $(date) ====="
