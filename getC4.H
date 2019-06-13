// Regularization filter and convection term
{
    // Calculate the C4 approximation for n+1 level
    volVectorField UProjF("UProjF", polyLaplaceFilter(UProj,deltaNf2,d1,d2) );
    UProjF.correctBoundaryConditions();

    // filtered flux
    surfaceScalarField phiF("phiF", polyLaplaceFilter(phiProj,deltaNf2,d1,d2));

    // Fist term: (C(us_f) uc_f)
    convection = explConvection(phiF,UProjF);

    // Add the second term Filt(C(us_f) uc')
    volVectorField Cint = explConvection(phiF, (UProj-UProjF));
    Cint.correctBoundaryConditions();

    convection += polyLaplaceFilter(Cint,deltaNf2,d1,d2);

    // Add the third term Filt(C(us') uc_f)
    Cint = explConvection((phiProj-phiF),UProjF);
    Cint.correctBoundaryConditions();

    convection += polyLaplaceFilter(Cint,deltaNf2,d1,d2);
}