/*---------------------------------------------------------------------------*\
  Selective Frequency Damping (SFD) as an fvOption.

  Encapsulated formulation (Jordi, Cotter & Sherwin, IJNMF 2014).

  Adds to the momentum equation the forcing:

      F = -chi * (U - Ubar)

  where Ubar is a temporally low-pass-filtered copy of U that satisfies

      dUbar/dt = (U - Ubar) / Delta

  The combination drives U toward a steady "base flow" of NS, killing any
  oscillatory mode of angular frequency omega satisfying

      chi   > |omega| / 2
      Delta > 1 / |omega|

  Encapsulated time-update of the filter (analytic, unconditionally stable):

      Ubar^{n+1} = Ubar^n + (1 - exp(-dt/Delta)) * (U^{n+1} - Ubar^n)

  Dictionary entry (system/fvOptions):

      SFD1
      {
          type            selectiveFrequencyDamping;
          libs            ("libSFDOption.so");
          active          true;
          selectionMode   all;
          fields          (U);        // velocity field to filter (the only field knob)
          chi             0.5;        // gain         [1/s]
          Delta           10.0;       // filter width [s]
      }
\*---------------------------------------------------------------------------*/

#ifndef SFDOption_H
#define SFDOption_H

#include "cellSetOption.H"
#include "volFields.H"

namespace Foam
{
namespace fv
{

class SFDOption
:
    public fv::cellSetOption
{
    // Private data

        //- Name of velocity field to filter (default "U")
        word UName_;

        //- SFD gain coefficient chi [1/s]
        scalar chi_;

        //- SFD filter width Delta [s]
        scalar Delta_;

        //- Filtered (target base-flow) field, registered to mesh
        autoPtr<volVectorField> UbarPtr_;

        //- Last simulation time at which Ubar was updated; gate against
        //  multiple updates per timestep across PIMPLE outer iterations
        scalar lastUpdateTime_;


    // Private member functions

        //- Construct Ubar field: read from disk if present (restart),
        //  else initialise from U
        void initialiseUbar();

        //- No copy
        SFDOption(const SFDOption&) = delete;
        void operator=(const SFDOption&) = delete;


public:

    //- Runtime type information
    TypeName("selectiveFrequencyDamping");


    // Constructors

        SFDOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~SFDOption() = default;


    // Member Functions

        // Sources (incompressible — pimpleFoam, simpleFoam, pisoFoam)

            //- Add forcing -chi*(U - Ubar) to the momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // Sources (compressible)

            //- Compressible variant: same forcing, multiplied by rho
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // Correction

            //- Update Ubar via encapsulated analytic step.
            //  Called by host solver via fvOptions.correct(U).
            virtual void correct(volVectorField& field);


        // IO

            //- Re-read parameters
            virtual bool read(const dictionary& dict);
};


} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
