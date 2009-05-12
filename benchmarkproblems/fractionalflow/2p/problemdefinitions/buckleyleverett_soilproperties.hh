// $Id: test_fractionalflow_soilproperties.hh 1475 2009-03-18 14:05:01Z markus $
#ifndef BUCKLEYLEVERETT_SOILPROPERTIES_HH
#define BUCKLEYLEVERETT_SOILPROPERTIES_HH

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** soil class for the buckley-leverett problem */

template<class Grid, class Scalar>
class BuckleyLeverettSoil: public Matrix2p<Grid, Scalar>
{
public:
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld, numEq = 1
    };
typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Matrix2p<Grid, Scalar>::modelFlag modelFlag;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return constPermeability_;
    }
    virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.2;
    }

    virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.2;
    }

    virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.2;
    }

    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {

        std::vector<double> param(2);

        //linear parameters
        param[0] = 0.; // minimal capillary pressure
        param[1] = 0.; // maximal capillary pressure

        //Brooks-Corey parameters
        //              param[0] = 0.; // lambda
        //              param[1] = 0.; // entry-pressure

        return param;
    }

    virtual modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid, Scalar>::linear;
    }

    BuckleyLeverettSoil()
    :Matrix2p<Grid,Scalar>()
    {
        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-7;
        }
    }

private:
    FieldMatrix constPermeability_;

};

} // end namespace
#endif
