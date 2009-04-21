/// $Id: co2_soilproperties.hh 610 2008-09-18 17:04:34Z melanie $
#ifndef TIMEDEPPROBLEM_SOIL
#define TIMEDEPPROBLEM_SOIL

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class Grid, class Scalar>
class NimasSoil: public Matrix2p<Grid, Scalar>
{
public:
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype CoordScalar;

    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

    typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    const FieldMatrix<CoordScalar,dim,dim> &K (const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        //        if (x[0] < 0.035)
        return permeability_;
        //        else
        //            return permeability2_;
    }

    double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return 0.42;
    }

    double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
    {
        return 0.0238;
    }

    double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
    {
        return 0.01;
    }

    virtual double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return     (800 /* spec. heat cap. of sediment */
                    * 2650 /* density of sediment */
                    * (1-porosity(x, e, xi)));
    }

    virtual double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
    {
        static const double ldry = 0.32;
        static const double lsat = 2.7;
        double l_wurz, Sw;
        Sw = sat;
        if(Sw<0.0) Sw = 0.0; /* ACHTUNG Regularisierung */
        if(Sw>1.) Sw = 1.; /* ACHTUNG Regularisierung */

        l_wurz = ldry + sqrt(Sw)*(lsat-ldry);

        if(isnan(l_wurz)) {
            std::cout <<"isnan heatcondwurzel \n"<<std::endl;
            l_wurz = 0.0;
        }
        return(l_wurz);
    }

    std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
    {
        std::vector<double> param(5);
        const double vgN = 7.06;

        //Van Genuchten parameters
        param[0] = 1 - 1/vgN; // m
        param[1] = vgN; // n
        param[2] = 0.5; // epsilon
        param[3] = 1.0/3.0; // gamma
        param[4] = 0.001019; // alpha

        return param;
    }

    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return Matrix2p<Grid,Scalar>::van_genuchten;
    }

    NimasSoil()
        :Matrix2p<Grid,Scalar>(),
         permeability_(0.0),
         permeability2_(0.0)
    {
        for (int i = 0; i < dim; i++){
            permeability_[i][i] = 0.5e-12;
            permeability2_[i][i] = 1.0e-13;
        }
    }

private:
    Dune::FieldMatrix<Scalar,dim,dim> permeability_;
    Dune::FieldMatrix<Scalar,dim,dim> permeability2_;
};
} // end namespace
#endif
