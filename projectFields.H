//Project velocity and mass flux and pressure to (n+1)th time step -----------//
volVectorField UProj(
    "UProj",
    ((1+kappa)*U.oldTime() - kappa*U.oldTime().oldTime())
);
UProj.correctBoundaryConditions();


surfaceScalarField phiProj(
   "phiProj",
   ((1+kappa)*phi.oldTime() - kappa*phi.oldTime().oldTime())
);

volScalarField pProj(
    "pProj",
    ((1+kappa)*p.oldTime() - kappa*p.oldTime().oldTime())
);
pProj.correctBoundaryConditions();


// filter residual of projected pressure (reduce checker-boards)
// demonstrated additional marginal dissipation
//{
//    volScalarField pProjF(polyLaplaceFilter(pProj,deltaNf2,d1,d2));
//    volScalarField pDiffF(polyLaplaceFilter((pProj-pProjF),deltaNf2,d1,d2));
//    pProj = (pProjF + pDiffF);
//    pProj.correctBoundaryConditions();
//}
//----------------------------------------------------------------------------//
