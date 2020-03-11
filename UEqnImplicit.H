// Project velocity into next time step --------------------------------------//

MRF.correctBoundaryVelocity(U);

tmp<volVectorField> tUe
(
    new volVectorField
    (
        IOobject
        (
            "Ue",
            U.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((1+kappa)*U.oldTime() - kappa*U.oldTime().oldTime())
    )
);
volVectorField& Ue = tUe.ref();


tmp<surfaceScalarField> tPhie
(
    new surfaceScalarField
    (
        IOobject
        (
            "phie",
            phi.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((1+kappa)*phi.oldTime() - kappa*phi.oldTime().oldTime())
    )
);
surfaceScalarField& phie = tPhie.ref();


fvVectorMatrix UEqn
(
      fvm::Sp((kappa + 0.5)/dt, U)
    + MRF.DDt(U)
    - fvm::laplacian(turbulence->nuEff(), U, "laplacian(nuEff,U)")
    ==
    (2*kappa*U.oldTime() - (kappa-0.5)*U.oldTime().oldTime()) / dt
    + fvc::div(turbulence->nuEff()*dev(T(fvc::grad(Ue, "grad(U)"))),
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


if(regOn)
{
    // Because of the non-iterative nature of the algorithy;
    // it is important to enforce the continuity error of extrapolated phi
    // is smaller than corrected phi;

    // RegModel->calcContinuityError(phie);

    solve(
        UEqn
        ==
        - fvc::grad(pProj, "grad(p)")   // pressure gradient
        - (
              fvc::div(phie,Ue, "div(phi,U)")
            - Ue * fvc::div(phie)
            - RegModel->convectionResidual(phie, Ue)
          )
    );
}
else
{
    solve(
        UEqn
        ==
        - fvc::grad(pProj, "grad(p)")   // pressure gradient
        - (
              fvc::div(phie,Ue, "div(phi,U)")
            - Ue * fvc::div(phie)
          )
    );
}

// clear temporary extrapolated fields
tUe.clear();
tPhie.clear();

phi = fvc::interpolate(U) & mesh.Sf(); // update mass flux

// ---------------------------------------------------------------------------//
