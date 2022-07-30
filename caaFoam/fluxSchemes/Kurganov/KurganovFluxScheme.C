#include "KurganovFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
defineTypeNameAndDebug(KurganovFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, KurganovFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KurganovFluxScheme::KurganovFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
:
    fluxScheme(typeName, dict),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),

    pos_(surfaceScalarField
    (
        IOobject
        (
            "pos",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("pos", dimless, 1.0)
    )),


    neg_(surfaceScalarField
    (
        IOobject
        (
            "neg",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("neg", dimless, -1.0)
    ))

{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

KurganovFluxScheme::~KurganovFluxScheme()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::KurganovFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{


    // --- upwind interpolation of primitive fields on faces

        surfaceScalarField rho_pos
        (
            "rho_pos",
            fvc::interpolate(rho_, pos_, "reconstruct(rho)")
        );
        surfaceScalarField rho_neg
        (
            "rho_neg",
            fvc::interpolate(rho_, neg_, "reconstruct(rho)")
        );

        surfaceVectorField rhoU_pos
        (
            "rhoU_pos",
            fvc::interpolate(rhoU_, pos_, "reconstruct(U)")
        );
        surfaceVectorField rhoU_neg
        (
            "rhoU_neg",
            fvc::interpolate(rhoU_, neg_, "reconstruct(U)")
        );

	const volScalarField& psi = thermo_.psi();
        volScalarField rPsi(1.0/psi);
        surfaceScalarField rPsi_pos
        (
            "rPsi_pos",
            fvc::interpolate(rPsi, pos_, "reconstruct(T)")
        );
        surfaceScalarField rPsi_neg
        (
            "rPsi_neg",
            fvc::interpolate(rPsi, neg_, "reconstruct(T)")
        );

	volScalarField ee("ee",rhoE_/rho_-0.5*magSqr(U_));	
        surfaceScalarField e_pos
        (
            "e_pos",
            fvc::interpolate(ee, pos_, "reconstruct(T)")
        );
        surfaceScalarField e_neg
        (
            "e_neg",
            fvc::interpolate(ee, neg_, "reconstruct(T)")
        );


        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);    

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh_.Sf());

        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh_.Sf());

        volScalarField c(sqrt(thermo_.Cp()/thermo_.Cv()*rPsi));

        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            fvc::interpolate(c, pos_, "reconstruct(T)")*mesh_.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            fvc::interpolate(c, neg_, "reconstruct(T)")*mesh_.magSf()
        );

	dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );

       surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );




        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));


        surfaceScalarField aSf("aSf", am*a_pos);

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

//        surfaceScalarField Uf_("Uf_", a_pos*U_pos+a_neg*U_neg);	

//	this->save(facei, patchi, aOwn*UOwn + aNei*UNei, Uf_);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

	surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
	surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

	amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));


	phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;


	phiUp= (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
		+ (a_pos*p_pos + a_neg*p_neg)*mesh_.Sf();

	phiEp= aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
		+ aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
		+ aSf*p_pos - aSf*p_neg;


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

