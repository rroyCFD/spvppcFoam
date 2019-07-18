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

#include "regularizationModel.H"
#include "TaylorGreenVortex2D.H"
#include "TaylorGreenVortex3D.H"

#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H" // for input/ouput format control

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    pisoControl piso(mesh, "SPVPPC");

    Info << "Non-orthogonal correctors: " <<  piso.nNonOrthCorr() << endl;

    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    dimensionedScalar dt = runTime.deltaT();
    const scalar kappa(0.5); // 2nd order explicit Adam-Bashforth scheme

    regularizationModel C4Regularization(U, phi, pp, pRefCell, pRefValue);
    C4Regularization.setRegOn(regOn);

//    TaylorGreenVortex2D TGV2D(U, phi, p, pRefCell);
//    TGV2D.setPrecision(8); // set log data precision

//    // initialize  time=0 fields as analytical
//    TGV2D.setInitialFieldsAsAnalytical();

//    // copy the analytical fields
//    if(runTime.startTime().value() == 0)
//    {
//        Info << "Time " << runTime.startTime().value()
//             << ": making U and p analytical; phi is interpolated." << endl;

//        U   = TGV2D.getUAnalytical();
//        phi = TGV2D.getPhiAnalytical();
//        p   = TGV2D.getPAnalytical();
//    }

    TaylorGreenVortex3D TGV3D(U, phi, p);
    TGV3D.setPrecision(8);
    TGV3D.setInitialFieldsAsAnalytical(U, phi, p);

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


        // Taylor-Green vortex
//        // calculate error fields and error norms
//        TGV2D.calcError();

//        // calculate Global Ek and epsilon values
//        TGV2D.calcGlobalProperties();

//        // write to log files
//        TGV2D.write();

        // calculate Global Ek and epsilon values
        TGV3D.calcGlobalProperties();

        // write to log files
        TGV3D.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
