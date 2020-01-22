// Non-orthogonal pressure corrector loop ------------------------------------//
MRF.makeRelative(phi);

adjustPhi(phi, U, p);

// Update the pressure BCs to ensure flux consistency
volScalarField rAU(1.0/UEqn.A());
constrainPressure(p, U, phi, rAU, MRF);


// RhieChow flux correction
if(RhieChow)
{
    tmp<surfaceVectorField> gradpDiff_f
    (
        (  (fvc::interpolate(fvc::grad(pProj, "grad(p)")) & ed)
          - faceGradient.snGrad(pProj)
        ) * ed
    );

    phi += fvc::interpolate(rAU)* (gradpDiff_f() & mesh.Sf());
    gradpDiff_f.clear();
}

// Non-orthogonal pressure corrector loop
pp.primitiveFieldRef() = 0.0;

while (piso.correctNonOrthogonal())
{
    fvScalarMatrix ppEqn
    (
      fvm::laplacian(pp, "laplacian(p)") == (fvc::div(phi)*((kappa + 0.5)/dt) )
    );

    // no correction to reference cell pressure
    scalar ppRefValue = pRefValue - getRefCellValue(pProj, pRefCell);
    // Info << "ppRefValue: " << ppRefValue << endl;

    ppEqn.setReference(pRefCell, ppRefValue); //
    ppEqn.solve(mesh.solver(pp.select(piso.finalNonOrthogonalIter())));

    if (piso.finalNonOrthogonalIter())
    {
        phi -= (ppEqn.flux() * (dt/(kappa + 0.5)) );
    }
}
#include "continuityErrs.H"

// Update the collocated velocity and pressure
p = (pProj+pp);
p.correctBoundaryConditions();

//check whether pp at refCell is smaller than solver tolerence!
if(pRefOn)
{
    // Info << "p at ref. cell: " << p[pRefCell] << tab
    //      << "pp at ref. cell: " << pp[pRefCell] << endl;
    Info << "p at ref. cell: " << getRefCellValue(p, pRefCell) << nl
         << "pp at ref. cell: " << getRefCellValue(pp, pRefCell) << endl;
}

U -= (fvc::grad(pp, "grad(p)") * (dt/(kappa + 0.5)) );

//U.correctBoundaryConditions(); // updated in laminar transport model
fvOptions.correct(U);
// ---------------------------------------------------------------------------//
