#
# airPower — central inputs
#
# All major simulation parameters are defined here.
# Access: inp.<CaseName>.<variable>
# blockMeshDict-related inputs are kept in blockMeshDict.
#

const inp = (

    # ======================================================================
    # DirectFlatPlate
    # ======================================================================
    DirectFlatPlate = (
        # Inflow parameters
        Uinf   = 12.417664315415696,    # chordwise velocity [m/s]
        Winf   = 20.3,                  # spanwise velocity [m/s]
        xInlet = 0.046694057793992,     # inlet position (distance from virtual LE) [m]

        # Top-boundary pressure polynomial (Casacuberta et al, 2022)
        pa4 = 0.004709401639645,
        pa3 = 0.059408065736933,
        pa2 = 0.245222700832888,
        pa1 = 0.657504955493668,
        pa0 = 2.006271563747756,

        # Fluid properties
        freeStreamViscosity = 1.456610719354608e-5,   # [m^2/s]

        # Parallel
        nProcs = 8,

        # Output settings
        outputFormat = "csv",   # csv | binary
    ),

    # ======================================================================
    # TunnelToCurvedPlate
    # ======================================================================
    TunnelToCurvedPlate = (
        # (to be added)
        nProcs = 12,
    ),

)
