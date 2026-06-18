#
# airPower — central inputs
#
# All major simulation parameters are defined here.
# Access: inp.<CaseName>.<variable>
#

const inp = (

    # ======================================================================
    # wallModulation — shared bump SHAPE (independent of placement coord)
    # Per-case positions live in inp.DFP.wallModulation and inp.TTCP.airfoilLE.wallModulation
    # ======================================================================
    wallModulation = (
        enabled = true,            # true to activate the bump
        mode    = :single,         # :single (parameters below) or :multiple (read from file, future)
        shape   = :esn,            # :sigmoidal or :esn (epsilon-skewed normal)
        A       = 1E-03,           # [m] height (positive=protrusion, negative=depression)

        # Sigmoidal shape parameters
        p       = 10,              # front steepness exponent (≥3 for C2)
        q       = 3,               # back steepness exponent (≥3 for C2)

        # ESN shape parameters
        epsilon = 0.0,             # skewness (-1 to 1)
        R       = 50,              # base-to-height ratio [-]; R*A is the physical base width [m]
        yTol    = 1e-5,            # [m] tail truncation tolerance
    ),

    # ======================================================================
    # DFP — DirectFlatPlateModule
    # ======================================================================
    DFP = (
        # Domain geometry 
        domainLength = 0.398763, # [m]
        domainHeight = 0.020046, # [m]

        # Inflow parameters
        Uinf   = 15.10,                 # chordwise velocity [m/s]
        Winf   = 18.7379374199,         # spanwise velocity [m/s]
        xInlet = 0.04684993006,         # inlet position (distance from virtual LE) [m]

        # Top-boundary pressure polynomial (Casacuberta et al, 2022)
        pa4 = 0.002343802682877,
        pa3 = 0.037696164789990,
        pa2 = 0.175152760791389,
        pa1 = 0.530290605857447,
        pa0 = 1.857425373851459,

        # Fluid properties
        freeStreamViscosity = 1.472e-5,   # [m^2/s]

        # Grid resolution (multipliers on base cell counts)
        gridXfactor = 2,   # streamwise: base (1) is 616 cells (144+24+48+24+280+96)
        gridYfactor = 1,   # wall-normal: base (1) is 160 cells (120 bottom + 40 top)

        # Parallel
        nProcs = 8,

        # Wall modulation position — coordinates relative to domain origin (x=0 at inlet)
        # Physical distance from virtual LE = xInlet + x.
        # Shape parameters (A, R, shape, etc.) live in the top-level inp.wallModulation.
        wallModulation = (
            # ESN center
            xCenter = 0.08833,         # [m] bump center
            # Sigmoidal extents
            xStart  = 0.155,           # [m] start of bump
            xPeak   = 0.157,           # [m] peak location
            xEnd    = 0.165,           # [m] end of bump
        ),

        # Output settings
        outputFormat = "bin",           # csv | binary
        wallExtrapolation = true,       # add wall point (u=v=w=0, p extrapolated) to output
    ),

    # ======================================================================
    # TTCP — TunnelToCurvedPlateModule
    # ======================================================================
    TTCP = (

        # ── Shared flow/physics parameters ──────────────────────────────
        flow = (
            freeStreamVelocityStreamwise  = 51.2284922,  # [m/s]
            freeStreamVelocitySpanwise    = 39.0256530,  # [m/s]
            freeStreamViscosity = 1.5582546e-5,          # [m^2/s]
        ),

        # ── TunnelCase-specific ─────────────────────────────────────────
        tunnel = (
            # Geometry
            tunnelInletHalfHeight  = 3.315,   # [m] (tunnel inlet half height, symmetry applied)
            tunnelOutletHalfHeight = 3.315,   # [m] (tunnel outlet half height, symmetry applied)
            tunnelLength           = 7.800,   # [m]

            # Airfoil geometry
            airfoilFile = "NACA0012.dat",  # file in io/airfoilGeometryData/
            chord       = 0.580716,    # [m]
            alphaDeg    = 5.0000,   # [deg]
            xCenter     = 0.0,      # [m]
            yCenter     = 0.0,      # [m]

            # Grid resolution
            Nx = 800,
            Ny = 360,

            # Turbulence
            turbulenceIntensity = 0.0007, # [-] (turbulenceIntensity = 0.0003 -> I = 0.03 %)
            turbLengthScale     = 0.009,  # [m]

            # BSpline control points (upper-wall curvature, symmetry applied)
            xcp1 = -1.950, ycp1 = 0.900,  # [m]
            xcp2 = 0.0,    ycp2 = 0.900,  # [m]
            xcp3 = 1.950,  ycp3 = 0.900,  # [m]
        ),

        # ── AirfoilLECase-specific ──────────────────────────────────────
        airfoilLE = (
            # Domain extent (expressed as dimensionless x/c)
            xiArch           = 0.01,   # x/c boundary between arch (curved LE) and straight blocks
            xiSuctionOutlet  = 0.35,   # x/c at the upper (suction-side) domain outlet
            xiPressureOutlet = 0.15,   # x/c at the lower (pressure-side) domain outlet

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
            NyBL     = 400,
            gradBL   = 200.0,  # BL wall-normal expansion ratio (cell size grows toward outer edge)
            gradArch = 18.0,   # arch streamwise grading (clusters cells toward the LE)

            # Export / post-processing
            outputFormat = "binary",       # csv | binary
            exportMode   = "partial",   # full | partial
            xiInlet      = 0.05,        # chord fraction for inlet boundary
            xiOutlet     = 0.35,        # chord fraction for outlet boundary
            exportHeight = 0.020,       # wall-normal from surface [m]
            wallExtrapolation = true,   # add wall point (u=v=w=0, p extrapolated) to output

            # Aerodynamic forces — integrate static pressure and wall shear
            # stress over a wall arclength range defined by chord fractions
            # [xiStart, xiEnd]. Consumed by
            # PreProcessing/src/source/aeroForces.jl.
            aeroForces = (
                xiStart = 0.069212,
                xiEnd   = 0.124212,
            ),

            # Boundary-layer integral metrics (consumed by
            # PreProcessing/src/source/blMetrics.jl)
            blMetrics = (
                # method — one of:
                #   "vorticityIntegralTrapezoidal"  trapezoidal sum between cell centres
                #   "vorticityIntegralMidpoint"     FV midpoint rule (cell extents, wall slice included)
                #   "maxProfile"                    max(u_tan) in the column
                #   "fixedHeight"                   u_tan at the topmost cell (≈ exportHeight)
                #   "pressureBernoulli"             Bernoulli with p_e = p_wall (Prandtl)
                method = "vorticityIntegralMidpoint",
            ),

            # Suppress mapped-pressure noise in the BL on suction/pressure outlets:
            # within `factor` × Blasius δ from the wall, the prescribed p is set
            # constant equal to the value at the first face beyond that band.
            clampOutletPressureBL       = false,
            clampOutletPressureBLFactor = 2.0,

            # Wall modulation position on the airfoil upper surface (chord fraction).
            # Shape parameters (A, R, shape, etc.) live in the top-level inp.wallModulation.
            wallModulation = (
                # ESN center
                xiCenter = 0.10,       # x/c on chord; converted to arc-length s at evaluation
                # Sigmoidal extents
                xiStart  = 0.13,       # x/c
                xiPeak   = 0.15,       # x/c
                xiEnd    = 0.20,       # x/c
            ),
        ),

        # ── Parallel ────────────────────────────────────────────────────
        nProcs = 8,
    ),

    # ======================================================================
    # VAL — Validation / internal preprocessing (not OpenFOAM case inputs)
    # ======================================================================
    VAL = (
        externalToScaling = (
            flowDataFile = "M3J_Case_Q_25_AOA_3.017_Re_2171263.mat",  # in io/airfoilFlowData/
            percentCrop  = 5.0,           # start cropping at this % of chord
            Nfit         = 4,             # polynomial order
            fitLaw       = :logarithmic,  # :monomial or :logarithmic
            valUe        = false,         # true: overlay the .mat reference Ue on
                                          #       the airfoilLE BLQuantitiesCompare plot
        ),
        valPIV  = false,  # true: overlay experimental PIV reference data on plots
        Gen     = 0,      # validation generation (reads from Validation/Gen{N}/)
        Case    = 0,      # validation case (reads from Validation/Gen{N}/Experimental/Case{M}/)

        # Spectral integral-boundary-layer (IBL) validation solver — shared by
        # the DFP (profile overlay) and TTCP (δ99/δ* comparison). When valBL is
        # true the relevant module runs solver_IBL_spectral (via MATLAB). The
        # wall-normal grid is Chebyshev: blNcheb nodes over [0, H] with median
        # node at blYi; blnx streamwise marching points. The domain height H is
        # taken per case from its own geometry — DFP: domainHeight, TTCP:
        # airfoilLE exportHeight (height of the exported data box).
        valBL   = false,    # true to run the IBL solver (overlay/comparison)
        blNcheb = 120,     # IBL wall-normal Chebyshev node count
        blYi    = 0.002,   # IBL Chebyshev node median [m] (cluster in the BL)
        blnx    = 1200,     # IBL streamwise marching points
    ),

)
