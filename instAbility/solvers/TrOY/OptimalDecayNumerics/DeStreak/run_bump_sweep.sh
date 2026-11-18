#!/bin/bash
# DeStreak base-flow sweep — generate flat + bump base flows via DirectFlatPlate
# in the PLAYGROUND airPower (never touches canonical). Each run: set the bump in
# the playground inputs.jl (line-scoped perl edits), run DFP, copy midPlane.csv
# to bumps/bf_<id>.csv. All at the playground's current grid resolution.
set -u
PG=/Users/jcasacubertap/OpenFOAM/cases/variants/playground/testDestreaking/airPower
DEST=/Users/jcasacubertap/OpenFOAM/cases/canonical/airPower/instAbility/TrOY/OptimalDecayNumerics/DeStreak
INP=$PG/inputs.jl
MID=$PG/PreProcessing/modules/directFlatPlateModule/postProcessing/midPlane.csv
export AIRPOWER_MATLAB=/Applications/MATLAB_R2025b.app/bin/matlab GKSwstype=nul
JULIA=/Users/jcasacubertap/.juliaup/bin/julia

set_field () { # anchor_regex  key  value
  perl -pi -e "if(/$1/){ s/(\\b$2\\s*=\\s*)[^,\\s]+/\${1}$3/ }" "$INP"
}

# id  A  R  xCenter  enabled   (R=50 gentle bumps; endTime bounded to 1500)
CFG=(
  "flat  1E-3  50  0.16  false"
  "p16   1E-3  50  0.16  true"
  "d16  -1E-3  50  0.16  true"
)

echo "===== DeStreak DFP sweep start: $(date) ====="
for row in "${CFG[@]}"; do
  read -r id A R xc en <<< "$row"
  if [ -s "$DEST/bumps/bf_${id}.csv" ]; then
    echo ""; echo "----- [$id] already present ($(du -h "$DEST/bumps/bf_${id}.csv"|cut -f1)) — skipping -----"
    continue
  fi
  echo ""; echo "----- [$id] A=$A R=$R xCenter=$xc enabled=$en  $(date +%T) -----"
  set_field "true to activate"       enabled  "$en"
  set_field "positive=protrusion"    A        "$A"
  set_field "base-to-height ratio"   R        "$R"
  set_field "m\\]  bump height"      A        "$A"
  set_field "ESN bump center"        xCenter  "$xc"
  # echo the resulting config lines for the log
  grep -nE "true to activate|positive=protrusion|base-to-height ratio|bump height|ESN bump center" "$INP"

  t0=$(date +%s)
  ( cd "$PG" && caffeinate -i "$JULIA" run.jl PreProcessing DirectFlatPlate all ) \
      > "$DEST/dfp_${id}.log" 2>&1
  rc=$?
  t1=$(date +%s)
  echo "[$id] DFP exit=$rc  elapsed=$((t1-t0))s"
  if [ -f "$MID" ]; then
    cp "$MID" "$DEST/bumps/bf_${id}.csv"
    echo "[$id] saved -> bumps/bf_${id}.csv ($(du -h "$DEST/bumps/bf_${id}.csv" | cut -f1))"
  else
    echo "[$id] ERROR: midPlane.csv not found — skipping"
  fi
done
echo ""; echo "===== DeStreak DFP sweep done: $(date) ====="
ls -la "$DEST/bumps/"
