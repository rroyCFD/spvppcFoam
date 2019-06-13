// Energy-conservative generalization of non-linear convection term
if(regularizationOrder == "C2")
{
    #include "getC2.H"
}
else if(regularizationOrder == "C4")
{
    #include "getC4.H"
}
else if(regularizationOrder == "C6")
{
    #include "getC6.H"
}
else
{
   convection = explConvection(phiProj,UProj); // default: regular projection
}
//----------------------------------------------------------------------------//