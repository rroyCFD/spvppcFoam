// Non-orthogonal pressure corrector loop ------------------------------------//
MRF.makeRelative(phi);

adjustPhi(phi, U, p);

// Update the pressure BCs to ensure flux consistency
volScalarField rAU(1.0/UEqn.A());
constrainPressure(p, U, phi, rAU, MRF);

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    fvScalarMatrix ppEqn
    (
      fvm::laplacian(pp, "laplacian(p)") == (fvc::div(phi)*((kappa + 0.5)/dt) )
    );

    // no correction to reference cell pressure
    scalar ppRefValue = pRefValue - getRefCellValue(pProj, pRefCell);
//    Info << "ppRefValue: " << ppRefValue << endl;

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
//Info << "pp at refCell: " << getRefCellValue(pp, pRefCell) << endl;
Info << "p at refCell:  " << getRefCellValue(p, pRefCell) << endl;

// //check whether pp at refCell is smaller than solver tolerence!
// if(pRefOn)
// {
//     Info << "p[pRefCell]: " << p[pRefCell] << tab
//          << "pp[pRefCell]: " << pp[pRefCell] << endl;
// }


U -= (fvc::grad(pp, "grad(p)") * (dt/(kappa + 0.5)) );
//U.correctBoundaryConditions(); // updated in laminar transport model
U.correctBoundaryConditions();
fvOptions.correct(U);
// ---------------------------------------------------------------------------//
