// Regularized convection term (6th order approximation)
{
    // Calculate the convection approximation for n level
    volVectorField UProjF("UProjF", polyLaplaceFilter(UProj,deltaNf2,d1,d2) );
    UProjF.correctBoundaryConditions();

    // filtered flux
    surfaceScalarField phiProjF("phiProjF",
                                polyLaplaceFilter(phiProj,deltaNf2,d1,d2));

    volVectorField UPrime( "UPrime", (UProj-UProjF));
    surfaceScalarField phiPrime("phiPrime", (phiProj-phiProjF));

    // Fist term: (C(us_f) uc_f)
    convection = explConvection(phiProjF,UProjF);

    // Second term: (C(us_f) uc')
    convection += explConvection(phiProjF, UPrime);

    // Third term:(C(us') uc_f)
    convection += explConvection(phiPrime,UProjF);

    // Fourth & final term: Filter(C(us') uc')
    volVectorField Cint = explConvection(phiPrime,UPrime);
    Cint.correctBoundaryConditions();
    convection += polyLaplaceFilter(Cint,deltaNf2,d1,d2);
}