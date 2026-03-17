#
# airPower — central inputs
#
# All major simulation parameters are defined here.
# blockMeshDict-related inputs are kept in blockMeshDict.
#

# ==========================================================================
# DirectFlatPlate
# ==========================================================================

# --- Inflow parameters ---
const Uinf   = 12.417664315415696    # chordwise velocity [m/s]
const Winf   = 20.3                  # spanwise velocity [m/s]
const xInlet = 0.046694057793992     # inlet position (distance from virtual LE) [m]

# --- Top-boundary pressure polynomial (Casacuberta et al, 2022) ---
const pa4 = 0.004709401639645
const pa3 = 0.059408065736933
const pa2 = 0.245222700832888
const pa1 = 0.657504955493668
const pa0 = 2.006271563747756

# --- Fluid properties ---
const freeStreamViscosity = 1.456610719354608e-5   # [m^2/s]

# --- Output settings ---
const outputFormat = "csv"   # csv | binary

# ==========================================================================
# TunnelToCurvedPlate
# ==========================================================================

# (to be added)
