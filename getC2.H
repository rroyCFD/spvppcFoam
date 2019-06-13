// Regularization filter and convection term
{
    // Calculate the C4 approximation for n+1 level
    volVectorField UProjF("UProjF", polyLaplaceFilter(UProj,deltaNf2,d1,d2) );
    UProjF.correctBoundaryConditions();

    // filtered flux
    surfaceScalarField phiF("phiF", polyLaplaceFilter(phiProj,deltaNf2,d1,d2));

    // Fist term: Filt(C(us_f) uc_f)
    convection = polyLaplaceFilter(explConvection(phiF,UProjF),deltaNf2,d1,d2);
}