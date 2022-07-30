Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& p = thermo.p();
volScalarField& e = thermo.he();

if (e.name() != "e")
{
    FatalErrorIn("cea")
        << "Only energy type internalEnergy supported."
        << nl << exit(FatalError);
}

const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.internalField()) > 0.0)
{
    inviscid = false;
}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "rhoBoundaryTypes.H"
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho(),
    rhoBoundaryTypes
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);


volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);


surfaceScalarField phi("phi", mesh.Sf() & fvc::interpolate(rhoU));


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
Info << endl;


Info<< "Creating flux scheme\n" << endl;
autoPtr<Foam::fluxScheme> flux_;
flux_ =
(
    Foam::fluxScheme::New(mesh.schemesDict(), thermo, rho, U, rhoU, rhoE)
);

	Foam::fluxScheme& flux = flux_();



Info << endl;

rho.correctBoundaryConditions();



surfaceVectorField phiUp
(
 IOobject
 (
  "phiUp",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimArea*rhoU.dimensions()*U.dimensions()
 );



surfaceScalarField phiEp
(
 IOobject
 (
  "phiEp",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimArea*rhoE.dimensions()*U.dimensions()
 );


surfaceVectorField Up
(
 IOobject
 (
  "Up",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimArea*U.dimensions()
 );

