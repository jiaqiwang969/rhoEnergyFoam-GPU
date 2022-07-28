// Evaluate Ducors sensor
   volScalarField divU  = fvc::div(U) ;
   volScalarField divU2 = divU*divU ;
   volScalarField rotU2 = magSqr(fvc::curl(U)) ;
// Internal field
//   forAll(U,icell)
//   {
//    ducSensor[icell] = max(-divU[icell]/Foam::sqrt(divU2[icell] + rotU2[icell] + epsilon),0.) ;
//   }
    ducSensor = max(-divU/Foam::sqrt(divU2 + rotU2 + eps  ),0.) ;
