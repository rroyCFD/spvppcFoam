## Application : Symmetry Preserving Velocity Projection, Pressure Correction Foam (spvppcFoam)

### Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

### Description
Transient solver for incompressible, turbulent flow, using the Symmetry Preserving, Velocity Projection, Pressure Correction (SPVPPC) algorithm. Implemented using the following the publications:

        Trias, F. X., Lehmkuhl, O., Oliva, A., Pérez-Segarra, C. D., & Verstappen, R. W. C. P. (2014).
        Symmetry-preserving discretization of Navier–Stokes equations on collocated unstructured grids.
        Journal of Computational Physics, 258, 246–267. https://doi.org/10.1016/j.jcp.2013.10.031

    and
        Verstappen, R. (2008).
        On restraining the production of small scales of motion in a turbulent channel flow.
        Computers & Fluids, 37(7), 887–897. https://doi.org/10.1016/J.COMPFLUID.2007.01.013

Sub-models include:

* turbulence modelling, i.e. laminar, RAS or LES
* run-time selectable MRF and finite volume options such as mean velocity, etc.


### Disclaiimer:

This application is built based on [OpenFOAM version-6](https://openfoam.org/release/6/). Please read the _About OpenFOAM_ section to learn more on OpenFOAM.

The application is free to use. The author neither provide any warranty nor shall be liable for any damage incurred from this application.



#### About OpenFOAM

OpenFOAM is the leading free, open source software for computational fluid dynamics (CFD), owned by the OpenFOAM Foundation and distributed exclusively under the [General Public Licence (GPL)](http://www.gnu.org/copyleft/gpl.html). The GPL gives users the freedom to modify and redistribute the software and a guarantee of continued free use, within the terms of the licence. To learn more visit [https://openfoam.org/](https://openfoam.org/)
