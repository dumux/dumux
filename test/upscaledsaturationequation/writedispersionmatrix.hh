#ifndef WRITEDISPERSIONMATRIX_HH
#define WRITEDISPERSIONMATRIX_HH

#include <boost/format.hpp>
#include "dumux/io/importfromdgf_level.hh"
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>


namespace Dune
{
template<class SoilType>
void writeDispersionMatrix(SoilType& soil,std::string data = "data")//, std::string dxx = "Dxxdata", std::string dyy = "Dyydata",  std::string saturation = "saturationdata")
{
    enum{dim = 2};
    typedef Dune::Matrix<Dune::FieldVector<double,1> > DataType;
    //    DataType Dxx(1,1);
    //    DataType Dyy(1,1);
    DataType collectedData(1,1);


    typedef Dune::UGGrid<dim> GridType;
    typedef typename GridType::LevelGridView GridView;

    Dune::importFromDGF<GridView,DataType>(collectedData, data);

    //    Dune::importFromDGF<DataType,GridView>(soil.getDispersionSat(), saturation);
    //    Dune::importFromDGF<DataType,GridView>(Dxx, dxx);
    //    Dune::importFromDGF<DataType,GridView>(Dyy, dyy);

    //    int lSize = Dxx.N();
    //    int cSize = Dyy.M();

    int lSize =collectedData.N();
    int cSize = collectedData.M()/(11);

    soil.getDispersion().setSize(lSize,cSize);
    soil.getDispersion()=0;
    soil.getDispersionSat().setSize(lSize,cSize);
    soil.getDispersionSat()=0;
    soil.getM().setSize(lSize,cSize);
    soil.getM()=0;
    soil.getMSat().setSize(lSize,cSize);
    soil.getMSat()=0;

    for (int i=0;i<lSize;i++)
        {
            for (int j=0;j<cSize;j++)
                {
                    //            soil.getDispersion()[i][j][0][0]=Dxx[i][j];
                    //            soil.getDispersion()[i][j][1][1]=Dyy[i][j];
                    soil.getDispersion()[i][j][0][0]=collectedData[i][j];
                    soil.getDispersion()[i][j][1][1]=collectedData[i][j+cSize];
                    soil.getDispersionSat()[i][j]= collectedData[i][j+2*cSize];
                    soil.getM()[i][j][0]=collectedData[i][j+3*cSize];
                    soil.getM()[i][j][1]=collectedData[i][j+4*cSize];
                    soil.getM()[i][j][2]=collectedData[i][j+5*cSize];
                    soil.getM()[i][j][3]=collectedData[i][j+6*cSize];

                    soil.getMSat()[i][j][0]= collectedData[i][j+7*cSize];
                    soil.getMSat()[i][j][1]= collectedData[i][j+8*cSize];
                    soil.getMSat()[i][j][2]= collectedData[i][j+9*cSize];
                    soil.getMSat()[i][j][3]= collectedData[i][j+10*cSize];
                    //            std::cout<<soil.getDispersionSat()[i][j]<<std::endl;
                    //            std::cout<<soil.getDispersion()[i][j]<<std::endl;
                }
        }
    return;
}
}
#endif
