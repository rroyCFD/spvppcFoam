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
    // Approach 1
    // solve(UEqn == - fvc::grad(pProj, "grad(p)"));   // pressure gradient

    //Approach 2
    // solve(UEqn == -fvc::reconstruct(fvc::snGrad(p)* mesh.magSf()));

    // Approach 3
    solve(UEqn == -getLSGrad(LHS, p)); // use least square minimization of pGrad

    // fvOptions.correct(U); // source corrrection applied after pressure solve
}

phi = fvc::interpolate(U) & mesh.Sf(); // update mass flux

// ---------------------------------------------------------------------------//
