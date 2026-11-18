#
# airPower — secondary settings
#
# Entries set once and rarely touched: meshing details, export options, solver
# and plot fine-tuning. Same block layout as inputs.jl; the two are merged into
# one `inp` (PreProcessing/src/auxiliary/loadInputs.jl), so these are read as
# inp.<CaseName>.<variable> exactly like the inputs. A key lives in exactly one
# of the two files.
#

const settings = (

    # ======================================================================
    # wallModulation — shared bump SHAPE
    # ======================================================================
    wallModulation = (
        mode    = :single,         # :single (parameters in inputs.jl) or :multiple (read from file, future)

        # Sigmoidal shape parameters
        p       = 10,              # front steepness exponent (≥3 for C2)
        q       = 3,               # back steepness exponent (≥3 for C2)

        # ESN tail
        yTol    = 1e-5,            # [m] tail truncation tolerance
    ),

    # ======================================================================
    # DFP — DirectFlatPlateModule
    # ======================================================================
    DFP = (
        wallModulation = (
            # Sigmoidal extents (used only by the sigmoidal shape; not part of
            # the scaling spec)
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

        # ── TunnelCase-specific ─────────────────────────────────────────
        tunnel = (
            # Geometry
            tunnelInletHalfHeight  = 0.900,   # [m] (tunnel inlet half height, symmetry applied)
            tunnelOutletHalfHeight = 0.900,   # [m] (tunnel outlet half height, symmetry applied)
            tunnelLength           = 7.800,   # [m]

            # Grid resolution
            Nx = 800,
            Ny = 240,

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
            # Wall-normal band heights
            hBL  = 0.280,   # [m] height of boundary-layer block band (wall-normal)
            hFar = 0.280,   # [m] height of far-field block band (beyond BL)

            # Block segmentation along C-path
            nSegSuction  = 6,      # blocks from suction outlet to xiArch (upper side)
            nSegArchUp   = 2,      # blocks from xiArch to LE (upper side)
            nSegArchLo   = 2,      # blocks from LE to xiArch (lower side)
            nSegPressure = 2,      # blocks from xiArch to pressure outlet (lower side)
            cosineArch   = true,   # cosine clustering of stations near the LE

            # Grid grading
            gradBL   = 200.0,  # BL wall-normal expansion ratio (cell size grows toward outer edge)
            gradArch = 18.0,   # arch streamwise grading (clusters cells toward the LE)

            # Export / post-processing
            outputFormat = "binary",       # csv | binary
            exportMode   = "partial",   # full | partial
            xiInlet      = 0.02,        # chord fraction for inlet boundary
            xiOutlet     = 0.50,        # chord fraction for outlet boundary
            exportHeight = 0.020,       # wall-normal from surface [m]
            wallExtrapolation = true,   # add wall point (u=v=w=0, p extrapolated) to output

            # Aerodynamic forces — integrate static pressure and wall shear
            # stress over a wall arclength range defined by chord fractions
            # [xiStart, xiEnd]. Consumed by
            # PreProcessing/src/source/aeroForces.jl.
            aeroForces = (
                xiStart = 0.05,
                xiEnd   = 0.30,
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
                method = "vorticityIntegralTrapezoidal",
            ),

            # Suppress mapped-pressure noise in the BL on suction/pressure outlets:
            # within `factor` × Blasius δ from the wall, the prescribed p is set
            # constant equal to the value at the first face beyond that band.
            clampOutletPressureBL       = false,
            clampOutletPressureBLFactor = 2.0,

            wallModulation = (
                # Sigmoidal extents
                xiStart  = 0.13,       # x/c
                xiPeak   = 0.15,       # x/c
                xiEnd    = 0.20,       # x/c
            ),
        ),
    ),

    # ======================================================================
    # VAL — Validation / internal preprocessing (not OpenFOAM case inputs)
    # ======================================================================
    VAL = (
        externalToScaling = (
            flowDataFile = "M3J_Case_Q_26.50_AoA_-2.9950_Re_2290411.mat",  # in io/airfoilFlowData/
            percentCrop  = 5.0,           # start cropping at this % of chord
            Nfit         = 4,             # polynomial order
            fitLaw       = :logarithmic,  # :monomial or :logarithmic
            valUe        = true,         # true: overlay the .mat reference Ue on
                                          #       the airfoilLE BLQuantitiesCompare plot
        ),

        # Spectral integral-boundary-layer (IBL) validation solver — shared by
        # the DFP (profile overlay) and TTCP (δ99/δ* comparison). When valBL is
        # true the relevant module runs solver_IBL_spectral (via MATLAB). The
        # wall-normal grid is Chebyshev: blNcheb nodes over [0, H] with median
        # node at blYi; blnx streamwise marching points. The domain height H is
        # taken per case from its own geometry — DFP: domainHeight, TTCP:
        # airfoilLE exportHeight (height of the exported data box).
        valBL   = true,    # true to run the IBL solver (overlay/comparison)
        blNcheb = 120,     # IBL wall-normal Chebyshev node count
        blYi    = 0.002,   # IBL Chebyshev node median [m] (cluster in the BL)
        blnx    = 1200,     # IBL streamwise marching points
    ),

    # ======================================================================
    # PostProcessing — MATLAB base-flow / perturbation analysis
    # ======================================================================
    PostProcessing = (
        plotUProfiles  = false,              # true: also draw profiles_u.png (non-dimensional u base flow + |u~|).
                                             #   Off for now; the w-vs-PIV validation figures are unaffected.
        valXcZoom      = [5.0, 25.0],        # zoom figure: PIV-station x/c window [%] ([] -> every station in the domain)
        valShareX      = true,               # w-validation figures: one x-limit per ROW across all stations, so the
                                             #   streamwise growth is readable. false -> per-panel autoscaling, which
                                             #   resolves each station's shape but makes every panel look alike.
        valYTop        = [],                 # w-validation figures: wall-distance limit [delta_0], shared by both rows.
                                             #   [] -> auto (1.5 x the perturbation extent over the plotted stations)
        # Plot window — ONE definition, shared by every PostProcessing figure
        # (base flow, perturbation shapes, profiles, Reynolds–Orr). See
        # PostProcessing/src/plotting/plotWindow.m.
        plot = (
            bufferFrac = 0.85,     # x-axis end, as a fraction of the streamwise extent: where the
                                   #   DeHNSSo outflow buffer begins, so the damped tail is not shown.
                                   #   1 -> full domain.
            # y-axis max above the wall. Each load mode is plotted in the unit system it is
            # COMPUTED in, and keeps it in the axes and legends too: the stability analysis
            # runs non-dimensionally, so loadFields stays in delta_0, while loadBF carries
            # the dimensional OpenFOAM grid. One key each, units stated, no conversion
            # between them — plotWindow picks the right one from loadMode.
            yMaxFields = 8.0,     # loadFields: wall-normal max in n/delta_0 ; [] -> yWallFrac
            yMaxBF     = 0.005,    # loadBF    : wall-normal max in METRES     ; [] -> yWallFrac
            yWallFrac  = 0.30,     # near-wall fraction of the domain height used when the above is []
        ),
        ro = (
            loadAnalysis = false,      # false: compute + save + plot; true: load io/output bundle & re-plot (no recompute)
            climManual   = false,      # colour range of the production figures (reynoldsOrrProd/Decomp):
                                       #   false: automatic, symmetric +-(99th percentile of |P|) in the plot window;
                                       #   true : fixed to clim below
            clim         = [-1.0, 1.0],  # [lower, upper] colour limits, in the figures' units (P/u_inf^3), when climManual
        ),
    ),

)
