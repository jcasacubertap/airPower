#
# airPower — central inputs
#
# All major simulation parameters are defined here.
# Access: inp.<CaseName>.<variable>
#

const inp = (

    # ======================================================================
    # DFP — DirectFlatPlateModule
    # ======================================================================
    DFP = (
        # Domain geometry 
        domainLength = 0.650, # [m]
        domainHeight = 0.020, # [m]

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
    # TTCP — TunnelToCurvedPlateModule
    # ======================================================================
    TTCP = (

        # ── Shared flow/physics parameters ──────────────────────────────
        flow = (
            freeStreamVelocity  = 24.84840467482210,        # [m/s]
            freeStreamViscosity = 1.456610719354608e-5,     # [m^2/s]
        ),

        # ── TunnelCase-specific ─────────────────────────────────────────
        tunnel = (
            # Geometry
            tunnelInletHalfHeight  = 0.900,   # [m] (tunnel inlet half height, symmetry applied)
            tunnelOutletHalfHeight = 0.900,   # [m] (tunnel outlet half height, symmetry applied)
            tunnelLength           = 7.800,   # [m]

            # Airfoil geometry
            airfoilFile = "M3J.dat",  # file in InputOutput/AirfoilData/
            chord       = 0.900,    # [m]
            alphaDeg    = -3.0,     # [deg]
            xCenter     = 0.0,      # [m]
            yCenter     = 0.0,      # [m]

            # Grid resolution
            Nx = 210,
            Ny = 60,

            # Turbulence
            turbulenceIntensity = 0.0003, # [-] (turbulenceIntensity = 0.0003 -> I = 0.03 %)
            turbLengthScale     = 0.009,  # [m]

            # BSpline control points (upper-wall curvature, symmetry applied)
            xcp1 = -1.950, ycp1 = 0.900,  # [m]
            xcp2 = 0.0,    ycp2 = 0.900,  # [m]
            xcp3 = 1.950,  ycp3 = 0.900,  # [m]
        ),

        # ── AirfoilLECase-specific ──────────────────────────────────────
        airfoilLE = (
            # Domain extent (expressed as dimensionless x/c)
            xiArch           = 0.02,   # x/c boundary between arch (curved LE) and straight blocks
            xiSuctionOutlet  = 0.50,   # x/c at the upper (suction-side) domain outlet
            xiPressureOutlet = 0.08,   # x/c at the lower (pressure-side) domain outlet

            # Wall-normal band heights
            hBL  = 0.280,   # [m] height of boundary-layer block band (wall-normal)
            hFar = 0.280,   # [m] height of far-field block band (beyond BL)

            # Block segmentation along C-path
            nSegSuction  = 6,      # blocks from suction outlet to xiArch (upper side)
            nSegArchUp   = 2,      # blocks from xiArch to LE (upper side)
            nSegArchLo   = 2,      # blocks from LE to xiArch (lower side)
            nSegPressure = 2,      # blocks from xiArch to pressure outlet (lower side)
            cosineArch   = true,   # cosine clustering of stations near the LE

            # Grid resolution
            NxTotal  = 800,
            NyBL     = 600,
            gradBL   = 200.0,  # BL wall-normal expansion ratio (cell size grows toward outer edge)
            gradArch = 18.0,   # arch streamwise grading (clusters cells toward the LE)

            # Export / post-processing
            outputFormat = "csv",       # csv | binary
            exportMode   = "partial",   # full | partial
            xiInlet      = 0.05,        # chord fraction for inlet boundary
            xiOutlet     = 0.50,        # chord fraction for outlet boundary
            exportHeight = 0.020,       # wall-normal from surface [m]
        ),

        # ── Mapping (TunnelCase → AirfoilLECase) ───────────────────────
        mapping = (
            tunnelTime = "163",    # converged time step to sample from
            zProbe     = 0.005,    # [m] for sampling
        ),

        # ── Parallel ────────────────────────────────────────────────────
        nProcs = 12,
    ),

)
