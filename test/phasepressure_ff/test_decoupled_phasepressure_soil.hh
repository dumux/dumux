// $Id$
#ifndef TEST_DECOUPLED_PHASEPRESSURE_SOIL_HH
#define TEST_DECOUPLED_PHASEPRESSURE_SOIL_HH

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class DecoupledPPTestSoil: public Matrix2p<Grid, Scalar>
{
public:
    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
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

            //linear parameters
            param[0] = 0.;
            param[1] = 0;

            //Brooks-Corey parameters
            param[0] = 2; // lambda
            param[1] = 0.; // entry-pressure

        return param;
    }

    virtual modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid, Scalar>::brooks_corey;
    }

    DecoupledPPTestSoil()
        :Matrix2p<Grid,Scalar>(),constPermeability(0)
    {}

private:
    FieldMatrix constPermeability;

};

} // end namespace
#endif /*MATRIXPROPERTIES*/
