#
# airPower — central inputs
#
# The parameters you set per case or per run. Entries set once and rarely
# touched (meshing details, export options, solver/plot fine-tuning) live in
# settings.jl, with the same block layout. The two are merged into one `inp`
# (PreProcessing/src/auxiliary/loadInputs.jl); a key lives in exactly one file.
# Access: inp.<CaseName>.<variable>
#

const inputs = (

    # ======================================================================
    # wallModulation — shared bump SHAPE (independent of placement coord)
    # Per-case positions live in inp.DFP.wallModulation and inp.TTCP.airfoilLE.wallModulation
    # ======================================================================
    wallModulation = (
        enabled = false,            # true to activate the bump
        shape   = :esn,            # :sigmoidal or :esn (epsilon-skewed normal)
        A       = 1E-03,           # [m] height (positive=protrusion, negative=depression)

        # ESN shape parameters
        epsilon = 0.0,             # skewness (-1 to 1)
        R       = 50,              # base-to-height ratio [-]; R*A is the physical base width [m]
    ),

    # ======================================================================
    # DFP — DirectFlatPlateModule
    # ======================================================================
    DFP = (
        # Domain geometry
        domainLength = 0.398763,    #0.650, # [m]
        domainHeight = 0.020046,     #0.020046, # [m]

        # Inflow parameters
        Uinf   = 15.10,                 # chordwise velocity [m/s]
        Winf   = 18.7379374199,                              # spanwise velocity [m/s]
        xInlet = 0.04684993006,                  # inlet position (distance from virtual LE) [m]

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

        # Wall modulation — coordinates relative to domain origin (x=0 at inlet).
        # Physical distance from virtual LE = xInlet + x. Shape parameters (R,
        # shape, esn/sigmoidal, etc.) live in the top-level inp.wallModulation.
        wallModulation = (
            # ── Bump specification: give EXACTLY TWO of {A, xCenter, Rek,
            #    AoverDstar}; the other two are solved at :prep from the IBL
            #    baseline (wallModulationScaling.jl). {A, xCenter} = direct
            #    geometry (no solve). Any Rek/AoverDstar target triggers the IBL
            #    inversion. A resolved A overrides the shared wallModulation.A.
            #    (ESN/xCenter only for now; sigmoidal scaling is future work.)
            A          = 0.001,        # [m]  bump height       (number, or nothing → solve)
            xCenter    = 0.08833,        # [m]  ESN bump center   (number, or nothing → solve)
            Rek        = nothing,           # [-]  roughness Re = u_b(y=A)·A/ν  (target, or nothing)
            AoverDstar = nothing,         # [-]  A / δ*(xCenter)              (target, or nothing)
        ),
    ),

    # ======================================================================
    # TTCP — TunnelToCurvedPlateModule
    # ======================================================================
    TTCP = (

        # ── Shared flow/physics parameters ──────────────────────────────
        flow = (
            freeStreamVelocityStreamwise  = 18.738329701443508,  # [m/s]
            freeStreamVelocitySpanwise    = 18.738329701443508,  # [m/s]
            freeStreamViscosity = 1.471748410882160e-05,          # [m^2/s]
        ),

        # ── TunnelCase-specific ─────────────────────────────────────────
        tunnel = (
            # Airfoil geometry
            airfoilFile = "M3J.dat",  # file in io/airfoilGeometryData/
            chord       = 0.900,    # [m]
            alphaDeg    = -2.995,   # [deg]
            xCenter     = 0.0,      # [m]
            yCenter     = 0.0,      # [m]
        ),

        # ── AirfoilLECase-specific ──────────────────────────────────────
        airfoilLE = (
            # Domain extent (expressed as dimensionless x/c)
            xiArch           = 0.02,   # x/c boundary between arch (curved LE) and straight blocks
            xiSuctionOutlet  = 0.50,   # x/c at the upper (suction-side) domain outlet
            xiPressureOutlet = 0.08,   # x/c at the lower (pressure-side) domain outlet

            # Grid resolution
            NxTotal  = 800,
            NyBL     = 400,

            # Wall modulation position on the airfoil upper surface (chord fraction).
            # Shape parameters (A, R, shape, etc.) live in the top-level inp.wallModulation.
            wallModulation = (
                # ESN center
                xiCenter = 0.15,       # x/c on chord; converted to arc-length s at evaluation
            ),
        ),

        # ── Parallel ────────────────────────────────────────────────────
        nProcs = 8,
    ),

    # ======================================================================
    # VAL — Validation / internal preprocessing (not OpenFOAM case inputs)
    # ======================================================================
    VAL = (
        valPIV  = false,  # true: overlay experimental PIV reference data on plots
        Gen     = 0,      # validation generation (reads from Validation/Gen{N}/)
        Case    = 0,      # validation case (reads from Validation/Gen{N}/Experimental/Case{M}/)
    ),

    # ======================================================================
    # PostProcessing — MATLAB base-flow / perturbation analysis
    # The dispatcher (run.jl) exports this block to a generated MATLAB config
    # (PostProcessing/inputs_gen.m) that run.m loads. `task` is normally set on
    # the command line — julia run.jl PostProcessing <task> — overriding the
    # default here.
    # ======================================================================
    PostProcessing = (
        task           = "importData",         # 'importData' | 'reynoldsOrrProdTerms'
        loadMode       = "loadBF",             # 'loadBF' | 'loadFields'
        caseType       = "TTCP",                # physical case: 'DFP' | 'TTCP' (base-flow source + validation-station mapping)
        fieldsFile     = "StabRes_cfi_flat_plate_Gen0_Case0.mat",      # (loadFields) from `instAbility DeHNSSo run`
        modeIdx        = Int[2,2],           # which spanwise modes to PLOT ([] -> all); the analysis always computes all.
                                             #   Repeated entries are de-duplicated on load (importData), so a mode is
                                             #   never drawn twice or listed twice in a legend.
        validation     = false,                # w-profile figures are always drawn; this picks how:
                                              #   false: stations in S*/delta_0, DeHNSSo alone (profiles_w_broad.png, amplitude_w.png)
                                              #   true : stations in x/c vs PIV (Gen/Case from the VAL block), needs the airfoil:
                                              #     profiles_w_validation.png       — zoom: stations ARE PIV stations, inside valXcZoom
                                              #     profiles_w_validation_broad.png — broad: stations across the whole numerical
                                              #       domain, with the PIV overlay wherever a station is close to a PIV one
    ),

)
