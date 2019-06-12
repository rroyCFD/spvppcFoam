// Project velocity into next time step --------------------------------------//

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
      fvm::Sp((kappa + 0.5)/dt, U)
    + MRF.DDt(U)
    - fvm::laplacian(turbulence->nuEff(), U, "laplacian(nuEff,U)")
    ==
    (2*kappa*U.oldTime() - (kappa-0.5)*U.oldTime().oldTime()) / dt
    - convection   // explicit convection term
    + fvc::div(turbulence->nuEff()*dev(T(fvc::grad(UProj, "grad(U)"))),
               "div(nuEff,devTGradU)") // deviatoric part
    + fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == - fvc::grad(pProj, "grad(p)"));   // pressure gradient)

    // fvOptions.correct(U); // source corrrection applied after pressure solve
}

phi = fvc::interpolate(U) & mesh.Sf(); // update mass flux

// ---------------------------------------------------------------------------//
