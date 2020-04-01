// Project velocity into next time step --------------------------------------//

MRF.correctBoundaryVelocity(U);

// Update U and phi with their respective extrapolated value
phi = ((1+kappa) * phi.oldTime() - kappa * phi.oldTime().oldTime());

U = ((1+kappa) * U.oldTime() - kappa * U.oldTime().oldTime());
// U.correctBoundaryConditions(); // updated at UEqn solve



fvVectorMatrix UEqn
(
      fvm::Sp((kappa + 0.5)/dt, U)
    + MRF.DDt(U)
    - fvm::laplacian(turbulence->nuEff(), U)
    ==
    (2 * kappa * U.oldTime() - (kappa-0.5) * U.oldTime().oldTime()) / dt
    - (
          fvc::div(phi,U)
        - U * fvc::div(phi)
      )
    + fvc::div(turbulence->nuEff() * dev(T(fvc::grad(U))))
    + fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

p = ((1+kappa) * p.oldTime() - kappa * p.oldTime().oldTime());
p.correctBoundaryConditions();

if(regOn)
{
    // Because of the non-iterative nature of the algorithm;
    // it is important to enforce the continuity error of extrapolated phi
    // is smaller than corrected phi;

    // RegModel->calcContinuityError(phi);

    solve(
        UEqn
        ==
        - fvc::grad(p)
        - (
            - RegModel->convectionResidual(phi, U)
          )
    );
}
else
{
    solve(
        UEqn
        ==
        - fvc::grad(p)

    );
}

// ---------------------------------------------------------------------------//
