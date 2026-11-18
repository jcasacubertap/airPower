#!/bin/bash
# DeStreak — overnight autonomous pipeline:
#   Phase 1  DirectFlatPlate base-flow sweep (flat + 4 protrusions + 1 depression)
#   Phase 2  gridgen each base flow (clean MATLAB session, common Malik extent)
#   Phase 3  corrected base-flow-oriented inverse design (Casacuberta et al. 2024)
#   Phase 4  summary figure
# Runs the PLAYGROUND airPower (private /tmp/DFP_*_pg work dirs — never touches canonical).
set -u
DEST=/Users/jcasacubertap/OpenFOAM/cases/canonical/airPower/instAbility/TrOY/OptimalDecayNumerics/DeStreak
cd "$DEST"
MB=/Applications/MATLAB_R2025b.app/bin/matlab
IDS="{'flat','p11','p15','p19','p23','d15'}"
BASIS="{'p11','p15','p19','p23'}"
stamp(){ echo "[$(date +%H:%M:%S)] $*"; }

stamp "PHASE 1: DFP base-flow sweep"
./run_bump_sweep.sh > sweep_dfp.log 2>&1
ncsv=$(ls bumps/bf_*.csv 2>/dev/null | wc -l | tr -d ' ')
stamp "PHASE 1 done — $ncsv/6 base flows in bumps/"
if [ "$ncsv" -lt 6 ]; then stamp "WARNING: fewer than 6 base flows; continuing with what exists"; fi

stamp "PHASE 2: gridgen sweep (60x1200, clean session)"
"$MB" -batch "destreak_gengrids_sweep($IDS,60,1200)" > gengrids_sweep.log 2>&1
ngrid=$(ls gsweep_*_60x1200.mat 2>/dev/null | wc -l | tr -d ' ')
stamp "PHASE 2 done — $ngrid grids"

stamp "PHASE 3: corrected base-oriented design (60x1200)"
"$MB" -batch "destreak_design_base($BASIS,60,1200);" > design_base.log 2>&1
stamp "PHASE 3 done"
grep -E "reach|c\*|sum I2|progress|per-bump" design_base.log | sed 's/^ *//'

stamp "PHASE 4: figure"
"$MB" -batch "destreak_plot_design(60,1200);" > plot_design.log 2>&1
stamp "ALL DONE"
