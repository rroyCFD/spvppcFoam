// Project velocity into next time step --------------------------------------//

MRF.correctBoundaryVelocity(U);

C4Regularization.update(); // update the convection term

fvVectorMatrix UEqn
(
      fvm::Sp((kappa + 0.5)/dt, U)
    + MRF.DDt(U)
    - fvm::laplacian(turbulence->nuEff(), U, "laplacian(nuEff,U)")
    ==
    (2*kappa*U.oldTime() - (kappa-0.5)*U.oldTime().oldTime()) / dt
    - C4Regularization.convectionTerm()   // explicit convection term
    + fvc::div(turbulence->nuEff()*dev(T(fvc::grad(C4Regularization.Ue(), "grad(U)"))),
               "div(nuEff,devTGradU)") // deviatoric part
    + fvOptions(U)
);


UEqn.relax();

fvOptions.constrain(UEqn);

volScalarField pProj
(
    "pProj",
    ((1+kappa)*p.oldTime() - kappa*p.oldTime().oldTime())
);
pProj.correctBoundaryConditions();

if (piso.momentumPredictor())
{
     solve(UEqn == - fvc::grad(pProj, "grad(p)"));   // pressure gradient
}

phi = fvc::interpolate(U) & mesh.Sf(); // update mass flux

// ---------------------------------------------------------------------------//
