/*---------------------------------------------------------------------------*\

Class
    Foam::KurganovFluxScheme

Description
    KurganovFluxScheme flux is the basic flux for rhoCentralFoam

SourceFiles
    KurganovFluxScheme.C


\*---------------------------------------------------------------------------*/
#ifndef KurganovFluxScheme_H
#define KurganovFluxScheme_H

#include "fluxScheme.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class KurganovFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class KurganovFluxScheme
:
    public fluxScheme
{

protected:
    // Protected data
        const fvMesh& mesh_;
        const psiThermo& thermo_;
        const volScalarField& rho_;
        const volVectorField& U_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;

        surfaceScalarField pos_;
        surfaceScalarField neg_;
public:
    //- Runtime type information
    TypeName("Kurganov");
        // Constructors

        //- Construct from components
        KurganovFluxScheme
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& rhoU,
            const volScalarField& rhoE
        );

    // Destructor

        ~KurganovFluxScheme();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up);

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //







