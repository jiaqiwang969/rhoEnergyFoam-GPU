#include "ausmPlusUpFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ausmPlusUpFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, ausmPlusUpFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ausmPlusUpFluxScheme::ausmPlusUpFluxScheme
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

ausmPlusUpFluxScheme::~ausmPlusUpFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ausmPlusUpFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{
    const volScalarField& p = thermo_.p();

    // Left and right states
    tmp< surfaceScalarField > phi_L(fvc::interpolate(U_, pos_, "reconstruct(U)") & mesh_.Sf());
    tmp< surfaceScalarField > phi_R(fvc::interpolate(U_, neg_, "reconstruct(U)") & mesh_.Sf());


//    // Flux relative to mesh movement
//    if (mesh_.moving())
//    {
//        fvc::makeRelative(phi_L.internalField(), U_);
//        fvc::makeRelative(phi_R.internalField(), U_);
//    }

    // Critical acoustic velocity (Liou 2006)
    tmp< volScalarField > gamma = thermo_.gamma();
    //tmp< volScalarField > H((rhoE_+p)/rho_);
    tmp< volScalarField > H
    (
//        (max(rhoE_/rho_,dimensionedScalar("0", rhoE_.dimensions()/rho_.dimensions(), 0.0)) +
//         max(p/rho_,dimensionedScalar("0", p.dimensions()/rho_.dimensions(), 0.0)))
        (max(rhoE_/rho_,dimensionedScalar("0", rhoE_.dimensions()/rho_.dimensions(), SMALL)) +
         max(p/rho_,dimensionedScalar("0", p.dimensions()/rho_.dimensions(), SMALL)))
    );

//    if (mesh_.moving())
//    {
//        H.internalField() -= 0.5*(U_&U_);
//        volVectorField Urel(U_);
//        Urel -= fvc::reconstruct(fvc::meshPhi(U_));
//        H.internalField() += 0.5*(Urel&Urel);
//    }
//    H().rename("H");
//    bound(H(), dimensionedScalar("H0", H().dimensions(), 0.0));

    tmp< volScalarField > c = sqrt(2.0*(gamma()-1.0)/(gamma()+1.0)*H());
    c->rename("c");
    gamma.clear();

    tmp< surfaceScalarField > c_L(fvc::interpolate(c(), pos_, "reconstruct(T)"));
    tmp< surfaceScalarField > c_R(fvc::interpolate(c(), neg_, "reconstruct(T)"));
    c->clear();
    c_L = c_L()*c_L()/max(c_L(), phi_L()/mesh_.magSf());
    c_R = c_R()*c_R()/max(c_R(),-phi_R()/mesh_.magSf());
    tmp< surfaceScalarField > c_face(min(c_L(),c_R()));
    c_L.clear();
    c_R.clear();

    // Critical Mach number
    tmp< surfaceScalarField > Mach_L(phi_L()/(c_face()*mesh_.magSf()));
    tmp< surfaceScalarField > Mach_R(phi_R()/(c_face()*mesh_.magSf()));

    // Coefficient to replace if statements
    tmp< surfaceScalarField > wM_l = pos(1.0-mag(Mach_L()));                    // |Ml| < 1 ? 1 : 0

    // Split Mach numbers
    tmp< surfaceScalarField > M_L1  =  0.5*(Mach_L()+mag(Mach_L()));
    tmp< surfaceScalarField > M_L2p =  0.25*sqr(Mach_L()+1.0);
    tmp< surfaceScalarField > M_L2m = -0.25*sqr(Mach_L()-1.0);

    // Mach number flux
    tmp< surfaceScalarField > Mach_plus_L =
    //      wM_l()       * M_L2p()                                              // beta = 0
          wM_l()         * M_L2p() * (1.0 - 2.0*M_L2m())                        // beta = 1/8
        + (1.0-wM_l()) * M_L1();

    // Pressure flux
    tmp< surfaceScalarField > p_plus_L =
           (1.0-wM_l())*M_L1()/(Mach_L()+VSMALL)
         + wM_l()      *(M_L2p()*(( 2.0 - Mach_L()) - 3.0*Mach_L()*M_L2m()));   //alpha = 3/16

    Mach_L.clear();

    wM_l.clear();
    M_L1.clear();
    M_L2p.clear();
    M_L2m.clear();

    // Coefficient to replace if statements
    tmp< surfaceScalarField > wM_r = pos(1.0-mag(Mach_R()));                    // |Mr| < 1 ? 1 : 0

    tmp< surfaceScalarField > M_R1   =  0.5*(Mach_R()-mag(Mach_R()));
    tmp< surfaceScalarField > M_R2m  = -0.25*sqr(Mach_R()-1.0);
    tmp< surfaceScalarField > M_R2p  =  0.25*sqr(Mach_R()+1.0);

    tmp< surfaceScalarField > Mach_minus_R =
    //      wM_r()         * M_R2m()                                            // beta = 0
          wM_r()         * M_R2m() * (1.0 + 2.0*M_R2p())                        // beta = 1/8
        + (1.0 - wM_r()) * M_R1();

    tmp< surfaceScalarField > p_minus_R =
           (1.0-wM_r())*M_R1()/(Mach_R()+VSMALL)
         + wM_r()      *(M_R2m()*((-2.0 - Mach_R()) + 3.0*Mach_R()*M_R2p()));   //alpha = 3/16

    Mach_R.clear();

    wM_r.clear();
    M_R1.clear();
    M_R2p.clear();
    M_R2m.clear();

    tmp< surfaceScalarField > p_L(fvc::interpolate(p, pos_, "reconstruct(rho)"));
    tmp< surfaceScalarField > p_R(fvc::interpolate(p, neg_, "reconstruct(rho)"));

    tmp< surfaceScalarField > rho_L = fvc::interpolate(rho_, pos_, "reconstruct(rho)");
    tmp< surfaceScalarField > rho_R = fvc::interpolate(rho_, neg_, "reconstruct(rho)");

    // Diffusive term
    tmp< surfaceScalarField > M_mean = 0.5*(phi_L()*phi_L() + phi_R()*phi_R())/(c_face()*mesh_.magSf()*c_face()*mesh_.magSf());   // Mean local Mach number
    tmp< surfaceScalarField > MDiff = -0.25*max((1.0-M_mean()),0.0)*(p_R() - p_L())/(0.5*(rho_L()+rho_R())*c_face()*c_face());    // 0 < K_p < 1, Liou suggest 0.25
    M_mean.clear();

    tmp< surfaceScalarField > Mach_1_2 = Mach_plus_L() + Mach_minus_R();
    Mach_1_2 =Mach_1_2 +  MDiff;
//    Mach_1_2.internalField() += MDiff();

    MDiff.clear();
    Mach_plus_L.clear();
    Mach_minus_R.clear();

    // Velocity difference or diffusive pressure
    // tmp< surfaceScalarField > pDiff = -0.25*p_plus_L()*p_minus_R()*(rho_L()+rho_R())*c_face()*((phi_R()-phi_L())/mesh_.magSf());  // 0 < Ku < 1, Liou suggest 0.75.

    phi_L.clear();
    phi_R.clear();
    rho_L.clear();
    rho_R.clear();

    surfaceScalarField p_1_2 = p_plus_L()*p_L() + p_minus_R()*p_R();
//    p_1_2 +=  pDiff()

    p_L.clear();
    p_R.clear();
    p_plus_L.clear();
    p_minus_R.clear();
//    pDiff.clear();

    tmp< surfaceVectorField > U_f = fvc::interpolate(U_, Mach_1_2(), "reconstruct(U)");

    surfaceScalarField rhoa_LR  = Mach_1_2()*c_face()*fvc::interpolate(rho_, Mach_1_2(), "reconstruct(rho)");
    surfaceVectorField rhoaU_LR = rhoa_LR*U_f();
    volScalarField ee("ee",rhoE_/rho_-0.5*magSqr(U_));
    surfaceScalarField rhoah_LR = rhoa_LR*(fvc::interpolate(ee, Mach_1_2(), "reconstruct(T)") + 0.5*magSqr(U_f())) + p_1_2*Mach_1_2()*c_face();
//    // NOTE: According to Liou enthalpy should be interpolated.
//    surfaceScalarField rhoah_LR = rhoa_LR*(fvc::interpolate(H(), Mach_1_2(), "reconstruct(T)"));

    // Face velocity for sigmaDotU (turbulence term)
    Up = U_f()*mesh_.magSf();
    U_f.clear();
    H.clear();
//    c_face.clear();

    phi = rhoa_LR*mesh_.magSf();
    phiUp = rhoaU_LR*mesh_.magSf() + p_1_2*mesh_.Sf();
    phiEp = rhoah_LR*mesh_.magSf();

/*
    if (mesh_.moving())
    {
        //phiEp += fvc::meshPhi(U_)*fvc::interpolate(p, Mach_1_2(), "reconstruct(T)");
        //phiEp += fvc::meshPhi(U_)*fvc::interpolate(rho_, Mach_1_2(), "reconstruct(rho)")*fvc::interpolate(p/rho_, Mach_1_2(), "reconstruct(T)"); //Ensure consistent interpolation with pressure term above
        //Not neccessary if using second last line
        phiEp += p_1_2 * fvc::meshPhi(U_); //Ensure consistent interpolation with pressure term above
    }
*/

/*
    forAll(mesh_.neighbour(),f)
    {
        if (Pstream::myProcNo() == 3 && mesh_.owner()[f] == 101106)
        {
            Pout<< "Face: " << f << " " << phi[f] << " " << rhoa_LR[f]*mesh_.magSf()[f] << " " << rhoaU_LR[f]*mesh_.magSf()[f] << " " << U_f()[f] << endl;
        }
        if (Pstream::myProcNo() == 3 && mesh_.neighbour()[f] == 101106)
        {
            Pout<< "Face: " << -f << " " << -phi[f] << " " << -rhoa_LR[f]*mesh_.magSf()[f] << " " << -rhoaU_LR[f]*mesh_.magSf()[f] << " " << U_f()[f] << endl;
        }
    }
*/

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

