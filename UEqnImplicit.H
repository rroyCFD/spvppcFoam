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


// C4Regularization.update(); // update the convection term

tmp<volVectorField> tconvTerm
(
    new volVectorField
    (
        IOobject
        (
            "convTerm",
            U.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("", dimAcceleration, vector::zero)
    )
);
volVectorField& convTerm = tconvTerm.ref();

if(regOn)
{
    convTerm = RegModel->getConvectionTerm(phie, Ue);
}
else
{
    convTerm = fvc::div(phie,Ue);
}


fvVectorMatrix UEqn
(
      fvm::Sp((kappa + 0.5)/dt, U)
    + MRF.DDt(U)
    - fvm::laplacian(turbulence->nuEff(), U, "laplacian(nuEff,U)")
    ==
    (2*kappa*U.oldTime() - (kappa-0.5)*U.oldTime().oldTime()) / dt
    - convTerm // dereferenced explicit convection term
    + fvc::div(turbulence->nuEff()*dev(T(fvc::grad(Ue, "grad(U)"))),
               "div(nuEff,devTGradU)") // deviatoric part
    + fvOptions(U)
);

// clear temporary extrapolated fields
tUe.clear();
tPhie.clear();
tconvTerm.clear();

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
