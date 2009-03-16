#ifndef WRITEDISPERSIONMATRIX_INTERFACE_HH
#define WRITEDISPERSIONMATRIX_INTERFACE_HH

#include <boost/format.hpp>
#include "dumux/io/importcorrection.hh"

namespace Dune
{
template<class CoarseParameterType, class Grid>
void writeCorrection(CoarseParameterType& parameters, std::string dataM = "dataM",
        std::string dataMSat = "dataMSat", std::string dataD = "dataDispersion",
        std::string dataDSat = "dataDispersionSat")
{
    enum
    {
        dim = Grid::dimension
    };
    typedef std::vector<std::vector<Dune::FieldVector<double,2*dim> > >
            DataTypeConv;
    typedef std::vector<std::vector<Dune::FieldVector<double,dim> > >
            DataTypeDisp;
    typedef std::vector<std::vector<Dune::FieldVector<double,1> > >
            DataTypeDispSat;

    DataTypeConv& getM = parameters.getFluxCorr().writeMatrix();
    DataTypeConv& getMSat = parameters.getFluxCorrSat().writeMatrix();
    DataTypeDisp& getDispersion = parameters.getDispersion().writeMatrix();
    DataTypeDispSat& getDispersionSat = parameters.getDispersionSat().writeMatrix();

    Dune::importFile<DataTypeConv,Grid>(getM, dataM);
    Dune::importFile<DataTypeConv,Grid>(getMSat, dataMSat);
    Dune::importFile<DataTypeDisp,Grid>(getDispersion, dataD);
    Dune::importFile<DataTypeDispSat,Grid>(getDispersionSat, dataDSat);

    return;
}
}
#endif
