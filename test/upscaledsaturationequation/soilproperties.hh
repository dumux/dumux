#ifndef MATRIXPROPERTIES
#define MATRIXPROPERTIES

#include <dumux/material/property_baseclasses.hh>
#include "dumux/material/randompermeability.hh"

namespace Dune
{

template<class Grid, class Scalar>
class HeterogeneousSoil: public Matrix2p<Grid, Scalar>
{
public:
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld, numEq = 1
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    virtual const FieldMatrix& K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
//        for(int i = 0; i < dim; i++)
//        {
//            constPermeability[i][i] = 1e-10;
//        }
//
//        return constPermeability;
        return randomPermeability.K(element);
    }

    virtual Scalar porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.2;
    }

    virtual Scalar Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const Scalar T = 283.15) const
    {
        return 0;
    }

    virtual Scalar Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const Scalar T = 283.15) const
    {
        return 0;
    }

    /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
     * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
    virtual Scalar heatCap(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 790 /* spec. heat cap. of granite */
        * 2700 /* density of granite */
        * porosity(globalPos, element, localPos);
    }

    virtual Scalar heatCond(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const Scalar sat) const
    {
        static const Scalar lWater = 0.6;
        static const Scalar lGranite = 2.8;
        Scalar poro = porosity(globalPos, element, localPos);
        Scalar lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
        Scalar ldry = pow(lGranite, (1-poro));
        return ldry + sqrt(sat) * (ldry - lsat);
    }

    virtual std::vector<Scalar> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const Scalar T = 283.15) const
    {
        std::vector<Scalar> param(2);

        //linear parameters
        param[0] = 0;
        param[1] = 0.;

        //Brooks-Corey parameters
//        param[0] = 2; // lambda
//        param[1] = 0.; // entry-pressure

        return param;
    }

    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
//        return Matrix2p<Grid,Scalar>::brooks_corey;
        return Matrix2p<Grid,Scalar>::linear;
    }

    HeterogeneousSoil(const Grid& grid,const char* name = "permeab.dat", const bool create = true)
    :randomPermeability(grid, grid.maxLevel(), name, create),constPermeability(0)
    {}

    ~HeterogeneousSoil()
    {}

private:
public:
    LevelRandomPermeability<Grid> randomPermeability;
    FieldMatrix constPermeability;
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
