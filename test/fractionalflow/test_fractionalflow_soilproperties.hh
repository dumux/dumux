// $Id$
#ifndef TEST_FRACTIONALFLOW_SOILPROPERTIES_HH
#define TEST_FRACTIONALFLOW_SOILPROPERTIES_HH

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class FractionalFlowTestSoil: public Matrix2p<Grid, Scalar>
{
public:
    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Matrix2p<Grid, Scalar>::modelFlag modelFlag;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        for(int i = 0; i < dim; i++)
        {
            constPermeability[i][i] = 1e-7;
        }
        return constPermeability;
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
//        if (globalPos[0]<=300)
//        {
            //linear parameters
            param[0] = 0.;
            param[1] = 100.;
//        }
//        else
//        {
            //Brooks-Corey parameters
//            param[0] = 0; // lambda
//            param[1] = 0.; // entry-pressure
//        }
        return param;
    }

    virtual modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        //        if (x[0]<=300)
        //            return 1;
        //        else
        return Matrix2p<Grid, Scalar>::linear;
    }

    FractionalFlowTestSoil()
        :Matrix2p<Grid,Scalar>(),constPermeability(0)//,permeability(g, name, create)
    {}

private:
    FieldMatrix constPermeability;

public:

    //    RandomPermeability<Grid> permeability;
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
