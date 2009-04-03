// $Id: test_fractionalflow_soilproperties.hh 1512 2009-03-31 09:26:10Z markus $
#ifndef TEST_FRACTIONALFLOW_SOILPROPERTIES_HH
#define TEST_FRACTIONALFLOW_SOILPROPERTIES_HH

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class InjectionSoil: public Matrix2p<Grid, Scalar>
{
public:
    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Matrix2p<Grid, Scalar>::modelFlag modelFlag;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos)
    {
        if (globalPos[0] >= innerLowerLeft_[0] && globalPos[0] <= innerUpperRight_[0] && globalPos[1]
            >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1])
            return innerPermeability_;
        else
            return outerPermeability_;
    }
    virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.2;
    }

    virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.;
    }

    virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.;
    }


    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {

        std::vector<double> param(2);

            //linear parameters
            param[0] = 0.;
            param[1] = 0.;

            //Brooks-Corey parameters
            param[0] = 2; // lambda
            param[1] = 10000.; // entry-pressure

        return param;
    }

    virtual modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid, Scalar>::brooks_corey;
    }

    InjectionSoil()
        :Matrix2p<Grid,Scalar>()
    {
        //define a lense
        innerLowerLeft_[0] = 30;
        innerLowerLeft_[1] = 30;
        innerUpperRight_[0] = 60;
        innerUpperRight_[1] = 45;

        outerPermeability_[0][0] = outerPermeability_[1][1] = 1e-10;
        outerPermeability_[0][1] = outerPermeability_[1][0] = 0;

        innerPermeability_[0][0] = innerPermeability_[1][1] = 1e-12;
        innerPermeability_[0][1] = innerPermeability_[1][0] = 0;
    }

private:
    FieldMatrix outerPermeability_;
    FieldMatrix innerPermeability_;
    GlobalPosition innerLowerLeft_;
    GlobalPosition innerUpperRight_;

};

} // end namespace
#endif /*MATRIXPROPERTIES*/
