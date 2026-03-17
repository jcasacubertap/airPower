#
# airPower — central inputs
#
# All major simulation parameters are defined here.
# Access: inp.<CaseName>.<variable>
# blockMeshDict-related inputs are kept in blockMeshDict.
#

const inp = (

    # ======================================================================
    # DFP — DirectFlatPlateModule
    # ======================================================================
    DFP = (
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
            freeStreamViscosity = 1.456610719354608e-5,      # [m^2/s]
            chord    = 900,      # [mm]
            alphaDeg = -3.0,     # [deg]
            xCenter  = 0.0,      # [mm]
            yCenter  = 0.0,      # [mm]
        ),

        # ── TunnelCase-specific ─────────────────────────────────────────
        tunnel = (
            # Geometry
            tunnelInletHeight  = 900.0,   # [mm] (tunnel half height)
            tunnelOutletHeight = 900.0,   # [mm]
            flatPlateLength    = 7800,    # [mm]
            z0 = -10,                     # [mm]
            z1 =  10,                     # [mm]
            seleAR            = 4,        # [-]
            seleHalfThickness = 10,       # [mm]
            lamTurbHeight     = 40,       # [mm]

            # Grid resolution
            Nx = 210,
            Ny = 60,

            # Turbulence
            turbulenceIntensity = 0.0003, # [-] (I = 0.03%)
            turbLengthScale     = 0.009,  # [m]

            # BSpline control points (upper-wall curvature)
            xcp1 = -1950, ycp1 = 900,
            xcp2 = 0,     ycp2 = 900,
            xcp3 = 1950,  ycp3 = 900,
        ),

        # ── AirfoilLECase-specific ──────────────────────────────────────
        airfoilLE = (
            # Domain extent (dimensionless x/c)
            xiArch           = 0.02,
            xiSuctionOutlet  = 0.50,
            xiPressureOutlet = 0.08,

            # Wall-normal band heights [mm]
            hBL  = 280,
            hFar = 280,

            # Spanwise [mm]
            zWidth = 2,

            # Block segmentation along C-path
            nSegSuction  = 6,
            nSegArchUp   = 2,
            nSegArchLo   = 2,
            nSegPressure = 2,
            cosineArch   = true,

            # Grid resolution
            nxTotal = 800,
            nyBL    = 600,
            nz      = 2,
            gradBL   = 200.0,
            gradArch = 18.0,

            # Export / post-processing
            outputFormat = "csv",       # csv | binary
            exportMode   = "partial",   # full | partial
            xiInlet      = 0.05,        # chord fraction for inlet boundary
            xiOutlet     = 0.50,        # chord fraction for outlet boundary
            exportHeight = 20,          # wall-normal from surface [mm]
        ),

        # ── Mapping (TunnelCase → AirfoilLECase) ───────────────────────
        mapping = (
            tunnelTime = "163",    # converged time step to sample from
            zProbe     = 0.005,    # z [m] for sampling
        ),

        # ── Parallel ────────────────────────────────────────────────────
        nProcs = 12,
    ),

)
