/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is developed to work with OpenFOAM. Same license protocol applies.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author:
    Rajib Roy
    University of Wyoming
    rroy@uwyo,edu, roy.rajib@live.com

Application
    spvppcFoam

Description
    Transient solver for incompressible, turbulent flow, using the Symmetry
    Preserving, Velocity Projection, Pressure Correction (SPVPPC) algorithm.
    Implemented using the following the papers:
        Trias, F. X., Lehmkuhl, O., Oliva, A., Pérez-Segarra, C. D., &
        Verstappen, R. W. C. P. (2014). Symmetry-preserving discretization of
        Navier–Stokes equations on collocated unstructured grids.
        Journal of Computational Physics, 258, 246–267.
        https://doi.org/10.1016/j.jcp.2013.10.031

    and
        Verstappen, R. (2008). On restraining the production of small scales
        of motion in a turbulent channel flow.
        Computers & Fluids, 37(7), 887–897.
        https://doi.org/10.1016/J.COMPFLUID.2007.01.013

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H" // for input/ouput format control

#include "orthogonalSnGrad.H"

#include "RegularizationModel.H"

#include "kineticEnergyAnalysis.H"

#include "TaylorGreenVortex.H"
#include "TaylorGreenVortex2D.H"
#include "TaylorGreenVortex3D.H"

#define watch(x) Info << (#x) << " is " << (x) << endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh, "SPVPPC");

    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // Rhie-Chow interpolation stuff
    const surfaceVectorField ed = mesh.delta()()/mag(mesh.delta()());
    Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

    dimensionedScalar dt = runTime.deltaT();
    const scalar kappa(0.5); // 2nd order explicit Adam-Bashforth scheme

    #include "createRegularization.H"
    #include "createGlobalProperties.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"
         dt = runTime.deltaT();

        // Velocity projection: Pressure correction steps
        {
            // Solve momentum conservation equation
            #include "UEqnImplicit.H"

            // Solve poisson pressure equation (mass conservation)
            #include "ppEqn.H"
        }

        // Solve for turbulence model related fields
        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        #include "writeGlobalProperties.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    #include "clearObjectPointers.H"

    return 0;
}
// ************************************************************************* //
