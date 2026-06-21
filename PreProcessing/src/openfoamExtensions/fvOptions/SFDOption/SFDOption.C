/*---------------------------------------------------------------------------*\
  Selective Frequency Damping (SFD) fvOption — implementation.
  See SFDOption.H for the method, formulation, and dictionary syntax.
\*---------------------------------------------------------------------------*/

#include "SFDOption.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(SFDOption, 0);
    addToRunTimeSelectionTable(option, SFDOption, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::SFDOption::initialiseUbar()
{
    IOobject UbarHeader
    (
        UName_ + "bar",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);

    if (UbarHeader.typeHeaderOk<volVectorField>(true))
    {
        Info<< "    SFD: reading existing field " << UName_ << "bar"
            << " (restart)" << endl;
        UbarPtr_.reset(new volVectorField(UbarHeader, mesh_));
    }
    else
    {
        Info<< "    SFD: initialising " << UName_ << "bar from "
            << UName_ << endl;
        UbarPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    UName_ + "bar",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                U
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::SFDOption::SFDOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    chi_(coeffs_.get<scalar>("chi")),
    Delta_(coeffs_.get<scalar>("Delta")),
    UbarPtr_(nullptr),
    lastUpdateTime_(-GREAT)
{
    // Velocity field to filter is taken from the standard `fields` entry
    // (single source of truth). SFD always filters the velocity vector;
    // there is no separate "U" knob.
    coeffs_.readEntry("fields", fieldNames_);
    UName_ = fieldNames_[0];
    fv::option::resetApplied();

    initialiseUbar();

    Info<< "    SFD parameters: chi = " << chi_
        << " [1/s], Delta = " << Delta_ << " [s]" << endl;

    if (chi_ <= 0)
    {
        FatalErrorInFunction
            << "SFD gain chi must be > 0 (got " << chi_ << ")"
            << exit(FatalError);
    }
    if (Delta_ <= 0)
    {
        FatalErrorInFunction
            << "SFD filter width Delta must be > 0 (got " << Delta_ << ")"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::SFDOption::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (!UbarPtr_)
    {
        return;
    }

    const volVectorField& Ubar = UbarPtr_();

    // Forcing on the physical RHS of the momentum equation:
    //
    //     F = -chi * (U - Ubar)
    //       = +chi * Ubar  -  chi * U
    //
    // Follows the canonical OpenFOAM implicit-feedback pattern (cf.
    // interRegionHeatTransfer: `eqn += htc*Tref - fvm::Sp(htc, T)`):
    //   (a) Explicit  +chi*Ubar -> added to the source vector
    //   (b) Implicit  -chi*U    -> subtract fvm::Sp(chi, U) so the term lands
    //                              as -chi*psi on the physical RHS (damping).
    //
    // NB: `eqn += fvm::Sp(chi, psi)` would put +chi*psi on the RHS, i.e.
    // anti-damping -> divergence (verified empirically: ||U-Ubar|| grows and
    // Cd goes negative). The sign below is the corrected one.

    // chi restricted to the selected cells (selectionMode all | cellSet | ...);
    // zero elsewhere, so SFD forces only the chosen region (e.g. a thin
    // near-wall band) and leaves the free-stream untouched.
    volScalarField::Internal chiField
    (
        IOobject
        (
            "sfd:chi", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    const labelList& cells = this->cells();
    forAll(cells, i)
    {
        chiField[cells[i]] = chi_;
    }

    eqn += chiField*Ubar() - fvm::Sp(chiField, eqn.psi());
}


void Foam::fv::SFDOption::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (!UbarPtr_)
    {
        return;
    }

    const volVectorField& Ubar = UbarPtr_();

    // rho*chi restricted to the selected cells (selectionMode); zero elsewhere.
    volScalarField::Internal rhoChi
    (
        IOobject
        (
            "sfd:rhoChi", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(rho.dimensions()/dimTime, Zero)
    );
    const labelList& cells = this->cells();
    forAll(cells, i)
    {
        rhoChi[cells[i]] = chi_*rho[cells[i]];
    }

    // Compressible: forcing per unit volume is rho*chi*(Ubar - U)
    eqn += rhoChi*Ubar() - fvm::Sp(rhoChi, eqn.psi());
}


void Foam::fv::SFDOption::correct(volVectorField& field)
{
    if (!UbarPtr_ || field.name() != UName_)
    {
        return;
    }

    // Gate: only update once per timestep, even if pimpleFoam calls
    // fvOptions.correct(U) inside multiple outer iterations.
    const scalar now = mesh_.time().value();
    if (now <= lastUpdateTime_ + SMALL)
    {
        return;
    }

    // Encapsulated filter update:  Ubar^{n+1} = Ubar^n + f*(U - Ubar^n),
    //   f = 1 - exp(-dt/Delta).
    // With local time stepping (ddt localEuler) the step varies per cell, so
    // use the reciprocal-local-timestep field "rDeltaT" (f_c uses dt_c =
    // 1/rDeltaT_c); otherwise fall back to the single global timestep.
    volVectorField& Ubar = UbarPtr_();
    vectorField& UbarI = Ubar.primitiveFieldRef();
    const vectorField& UI = field.primitiveField();

    if (mesh_.foundObject<volScalarField>("rDeltaT"))
    {
        const scalarField& rDeltaT =
            mesh_.lookupObject<volScalarField>("rDeltaT").primitiveField();
        forAll(UbarI, c)
        {
            const scalar f =
                scalar(1) - std::exp(-scalar(1)/(rDeltaT[c]*Delta_));
            UbarI[c] += f*(UI[c] - UbarI[c]);
        }
    }
    else
    {
        const scalar f =
            scalar(1) - std::exp(-mesh_.time().deltaTValue()/Delta_);
        UbarI += f*(UI - UbarI);
    }
    Ubar.correctBoundaryConditions();

    lastUpdateTime_ = now;

    // Diagnostic: ||U - Ubar||_2 over the SFD region (this->cells()).
    const labelList& cells = this->cells();
    const scalarField& V = mesh_.V();
    scalar num = 0, den = 0;
    forAll(cells, i)
    {
        const label c = cells[i];
        num += magSqr(UI[c] - UbarI[c])*V[c];
        den += V[c];
    }
    reduce(num, sumOp<scalar>());
    reduce(den, sumOp<scalar>());
    const scalar errNorm = (den > 0) ? std::sqrt(num/den) : scalar(0);

    Info<< "    SFD: ||U - Ubar||_2 = " << errNorm << endl;
}


bool Foam::fv::SFDOption::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("chi",   chi_);
        coeffs_.readEntry("Delta", Delta_);
        return true;
    }
    return false;
}


// ************************************************************************* //
