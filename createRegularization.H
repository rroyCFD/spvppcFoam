// RegularizationModel

RegularizationModel* RegModel;
if(regOn)
{
    RegModel = new RegularizationModel(
                            pp, pRefCell, pRefValue, piso.nNonOrthCorr());
    RegModel->setRegOn(regOn);
    RegModel->setFilterFieldsDivFree(filterFieldDivFree);
    Info << "Regularization Object Constructor" << endl;
}
else
{
    delete RegModel;
}

//------------------------------------------------------------------------//