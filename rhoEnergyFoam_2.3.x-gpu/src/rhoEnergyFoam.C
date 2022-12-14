////  =========                 |
//  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
//   \\    /   O peration     |
//    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
//     \\/     M anipulation  |
//-------------------------------------------------------------------------------
//License
//    This file is part of OpenFOAM.
//
//    OpenFOAM is free software: you can redistribute it and/or modify it
//    under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
//    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//    for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
//
//Application
//    rhoEnergyFoam 
//
//Description
//    Numerical solver for the solution of compressible shock-free flows.
//    The convective terms are discretized using Pirozzoli's scheme(JCP 2010),
//    density based splitting. The scheme allows conservation of the kinetic
//    energy in the inviscid incompressible limit.
//    midPoint interpolation must be selected in fvScheme dictionary. 
//    Viscous terms(Laplacians only) are evaluated directly, computing
//    the face normal gradients.
//    A third-order low-storage RK scheme is used for time integration.
//    The OpenFOAM labrary for turbulence models is included.
//      
//    Author: Davide Modesti (davide.modesti@uniroma1.it)
//    Last update 06/04/2017
//    Reference
//    D. Modesti and S. Pirozzoli A high-fidelity solver for turbulent compressible flows on unstructured meshes. Comput. & Fluids (2017)
//\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include <fstream>      // std::ofstream

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// AUSM useful functions
  surfaceScalarField m1 (surfaceScalarField rm,float sgn) // Return M_1
   {
   surfaceScalarField r1("r1", 0.5*(rm + sgn*Foam::mag(rm))) ; 
   return r1 ; 
   }
///
  surfaceScalarField m2 (surfaceScalarField rm,float sgn) // Return M_2
   {
   surfaceScalarField r2("r2", sgn*0.25*(rm + sgn)*(rm + sgn)) ;
   return r2 ; 
   }
///
  surfaceScalarField p5 (surfaceScalarField rm, float sgn,surfaceScalarField alpha) // Return p_5
   {
   surfaceScalarField r3("r3", Foam::neg(Foam::mag(rm)-1.0)*m2(rm,sgn)*( (sgn*2-rm) - sgn*16*alpha*rm*m2(rm,-sgn)) + Foam::pos(Foam::mag(rm)-1.0)*(1./rm*m1(rm,sgn))) ; 
   return r3; 
/*   
   if (Foam::mag(rm)<1.)
   {
    r = m2(rm,sgn)*( (sgn*2-rm) - sgn*16*alpha*rm*m2(rm,-sgn) ) ;
   }
   else
   {
    r = 1./rm*m1(rm,sgn) ;
   }
*/   
   }
///
  surfaceScalarField m4 (surfaceScalarField rm,float sgn ,float beta) // Return M_4
   {
   surfaceScalarField r4("r4", Foam::neg(Foam::mag(rm)-1.0)*m2(rm,sgn)*(1-sgn*16*beta*m2(rm,-sgn))+ Foam::pos(Foam::mag(rm)-1.0)*m1(rm,sgn)) ; 
   return r4 ; 
/*
   if (abs(rm)<1)
   {
    r = m2(rm,sgn)*( 1 - sgn*16*beta*m2(rm,-sgn) ) ;
   }
   else
   {
    r = m1(rm,sgn) ;
   }
*/   
   }
// Main
int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMeshNoClear.H"
    #include "createFields.H"
    #include "readThermophysicalProperties.H"
    #include "readTimeControls.H"
    #include "variables.H"
    //  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;
    Info<< "Start Timing = " << runTime.clockTimeIncrement() << " s"
        << nl << endl;
    while (runTime.loop()) //Start time loop
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;


//     Saving quantities at preavious time step
       rhoOld  = rho; 
       rhoUOld = rhoU; 
       rhoEOld = rhoE; 
//     RK Time step
       for (int cycle =0; cycle < rkCoeff.size(); cycle++)
       {
//     Speed of sound and Mach number
        c = Foam::sqrt(thermo.Cp()/thermo.Cv()/psi);
        Mach = U/c ;
//      Interpolated quantities at cell faces     
        surfaceScalarField rhoave = fvc:: interpolate(rho) ;
        surfaceVectorField Uave   = fvc:: interpolate(U)   ;
//      Flux at the intercell
        phi    = fvc:: interpolate(U)    & mesh.Sf() ;
        phit   = fvc:: interpolate(rhoU) & mesh.Sf() ; //flux in turbulent model
//      Enthalpy
        H    = (rhoE + p)/rho ;
//      Enthalpy at the intercell
        surfaceScalarField Have = fvc:: interpolate(H) ;
//      Pressure at the intercell
        surfaceScalarField pave = fvc::interpolate(p)  ;       
        #include "sensor.H" //Ducros sensor
//        if(pressArtDiff)
//        {
//         #include "AUSM.H"    // AUSM+up dissipation on pressure term, add dissipation on pave
//        }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AUSM.H start~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/        
//     AUSM+up dissipation on pressure term, add dissipation on pave
//     Add AUSM dissipation on pressure term
//     Reconstruct M_L and M_R at left and right cell states
//     M_L = fvc::interpolate(Mach, pos, "reconstruct(M)") & (mesh.Sf()/(mesh.magSf()));
//     M_R = fvc::interpolate(Mach, neg, "reconstruct(M)") & (mesh.Sf()/(mesh.magSf()));
//     Reconstruct p_L = p_j  p_R = p_j+1
       p_L = fvc::interpolate(p, pos, "reconstruct(p)");
       p_R = fvc::interpolate(p, neg, "reconstruct(p)");
//     Reconstruct rho, U, c
       rho_L= fvc::interpolate(rho, pos, "reconstruct(rho)") ;
       rho_R= fvc::interpolate(rho, neg, "reconstruct(rho)") ;
       U_L  = fvc::interpolate(U, pos, "reconstruct(U)")& (mesh.Sf()/(mesh.magSf())) ;
       U_R  = fvc::interpolate(U, neg, "reconstruct(U)")& (mesh.Sf()/(mesh.magSf())) ;
//     c_L  = fvc::interpolate(c, pos, "reconstruct(T)") ;
//     c_R  = fvc::interpolate(c, neg, "reconstruct(T)") ;
//     surfaceScalarField cave = fvc::interpolate(c) ;
       const   labelgpuList& own = mesh.owner();
       surfaceScalarField duc = fvc::interpolate(ducSensor) ;
       volScalarField gamma = thermo.Cp()/thermo.Cv();
       surfaceScalarField gamm = fvc::interpolate(gamma);
       volScalarField cts = Foam::sqrt(2.*(gamma-1.)/(gamma+1.)*H);
       surfaceScalarField css = fvc::interpolate(cts);
       surfaceScalarField cs_R=css*css/max(css,-U_R) ;
       surfaceScalarField cs_L=css*css/max(css, U_L) ;
       c12=fvc::interpolate(c);
//     c12=Foam::sqrt(c_L*c_R);//min(cs_L,cs_R) ;
//     c12=min(cs_L,cs_R) ;
       M_L  = U_L/c12;
       M_R  = U_R/c12;
       m2= (U_L*U_L+U_R*U_R)/(2.*c12*c12) ;
       m0= Foam::sqrt( min( pos, max(m2,minfty) ) );

       surfaceScalarField fa("fa", m0*(2.-m0));
       surfaceScalarField alpha("alpha", 3./16.*(-4.+5.*fa*fa));

       surfaceScalarField p5p("p5p",    p5 (M_L , 1.0 , alpha));
       surfaceScalarField p5m("p5m",    p5 (M_R ,-1.0 , alpha));

       surfaceScalarField dpr("dpr",    p5 (M_R , 1 , alpha) - p5m);
       surfaceScalarField dpl("dpl",    p5p  - p5(M_L, -1, alpha));

       surfaceScalarField pu("pu",     -ku*p5p*p5m*(rho_R + rho_L)*c12*fa*(U_R-U_L));
       surfaceScalarField dp12("dp12",  p_R*dpr - p_L*dpl);

       //Update p
       pave +=  Foam::pos(duc-ducLevelPress)*duc*(-0.5*(dp12) + pu);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AUSM.H end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/

//
//      Evaluate viscous terms
//
//      Divergence of the velocity        
        vecDivU.component(0) = fvc::interpolate(fvc::div(U));
        vecDivU.component(1) = fvc::interpolate(fvc::div(U));
        vecDivU.component(2) = fvc::interpolate(fvc::div(U));
//
// 
        volScalarField muEff(turbulence->muEff());
        surfaceScalarField  muave  = fvc::interpolate(muEff);//mu at cell faces
//
        volScalarField k("k", thermo.Cp()*muEff/Pr);//thermal diffusivity
//
        surfaceScalarField kave=fvc::interpolate(k);//k at cell faces. alphaEff=muEff/Prt
        //momentum viscous flux
        surfaceVectorField momVisFlux = muave*(fvc::snGrad(U)*mesh.magSf() -  2./3.*vecDivU*mesh.magSf());
        //energy viscous flux
        surfaceScalarField heatFlux =  kave*fvc::snGrad(T)*mesh.magSf();
        surfaceScalarField visWork  =  momVisFlux & Uave;
        enVisFlux = heatFlux + visWork ;
//
        // Total fluxes, Eulerian + viscous 
        surfaceScalarField rhoFlux   = rhoave*phi                                     ;
        momFlux                      = rhoave*Uave*phi + pave*mesh.Sf() - momVisFlux  ;
        enFlux                       = rhoave*Have*phi - enVisFlux                    ;
//
//        if(convArtDiff)
//        {
//         #include "AUSM_conv.H"    // AUSM+up dissipation on convective terms, add dissipation on rhoFlux,momFlux,enFlux
//        }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AUSM_conv.H start~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/
// AUSM+up dissipation on convective terms, add dissipation on rhoFlux,momFlux,enFlux
       volScalarField rhoH = rho*H ;
//
       surfaceVectorField rhoU_R=fvc::interpolate(rhoU,   neg, "reconstruct(U)") ;
       surfaceVectorField rhoU_L=fvc::interpolate(rhoU,   pos, "reconstruct(U)") ;
       surfaceScalarField rhoH_L=fvc::interpolate(rhoH, pos, "reconstruct(T)") ;
       surfaceScalarField rhoH_R=fvc::interpolate(rhoH, neg, "reconstruct(T)") ;
//
       scalar beta = 1./8. ;

       surfaceScalarField dml("dml",    m4(M_L,1.,beta) - m4(M_L,-1.,beta));
       surfaceScalarField dmr("dmr",    m4(M_R,1.,beta) - m4(M_R,-1.,beta));
       surfaceScalarField dp("dp",      -kp/fa*Foam::max(1-sigma*m2,0.0)*(p_R-p_L)/(rhoave*c12*c12));

       surfaceScalarField dm12("dm12",   dmr - dml - 2.*dp);
       surfaceScalarField m12("m12",     0.5*(M_L + M_R)  - 0.5*dm12);
       surfaceScalarField dw1("dw1",     (0.5*dm12 -Foam::mag(m12))*rho_L + (0.5*dm12 + Foam::mag(m12))*rho_R);
       surfaceScalarField dw2("dw2",     (0.5*dm12 -Foam::mag(m12))*rhoU_L.component(0) + (0.5*dm12 + Foam::mag(m12))*rhoU_R.component(0));
       surfaceScalarField dw3("dw3",     (0.5*dm12 -Foam::mag(m12))*rhoU_L.component(1) + (0.5*dm12 + Foam::mag(m12))*rhoU_R.component(1));
       surfaceScalarField dw4("dw4",     (0.5*dm12 -Foam::mag(m12))*rhoU_L.component(2) + (0.5*dm12 + Foam::mag(m12))*rhoU_R.component(2));
       surfaceScalarField dw5("dw5",     (0.5*dm12 -Foam::mag(m12))*rhoH_L + (0.5*dm12 + Foam::mag(m12))*rhoH_R);

       //Update convective fluxes
       rhoFlux -=  Foam::pos(duc-ducLevelConv)*0.5*c12*dw1*mesh.magSf();
       momFlux.component(0) = momFlux.component(0) - Foam::pos(duc-ducLevelConv)*0.5*c12*dw2*mesh.magSf();
       momFlux.component(1) = momFlux.component(1) - Foam::pos(duc-ducLevelConv)*0.5*c12*dw3*mesh.magSf();
       momFlux.component(2) = momFlux.component(2) - Foam::pos(duc-ducLevelConv)*0.5*c12*dw4*mesh.magSf();
       enFlux -= Foam::pos(duc-ducLevelConv)*0.5*c12*dw5*mesh.magSf();


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AUSM_conv.H end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/
// 
        volScalarField rhoFl = fvc::div(rhoFlux);
        volVectorField momFl = fvc::div(momFlux);
        volScalarField enFl  = fvc::div(enFlux);
        // RK sub-step
        rho  = rhoOld + rkCoeff[cycle]*runTime.deltaT()*(
                - rhoFl);         
//
        rhoU = rhoUOld + rkCoeff[cycle]*runTime.deltaT()*(
                - momFl);
//
        rhoE = rhoEOld + rkCoeff[cycle]*runTime.deltaT()*(
                -enFl);
        //Update primitive variables and boundary conditions
        U.dimensionedInternalField() = rhoU.dimensionedInternalField() / rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField()         = rho.boundaryField()*U.boundaryField();

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        //Thermodinamic library
        thermo.correct();
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
//
        p.dimensionedInternalField() = rho.dimensionedInternalField() / psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() = psi.boundaryField()*p.boundaryField(); //psi=1/(R*T)
        runTime.write();
        turbulence->correct(); //turbulence model
       }//end of RK time integration
//
//        #include "diagnostics.H" //print tke on diagnostics.dat
        #include "step.H"        //Evaluate Courant Number
        #include "setDeltaT.H"   //Andjust time step
//
// 
    }
        runTime.write();
    Info<< "Start Timing = " << runTime.clockTimeIncrement() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
