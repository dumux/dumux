/// $Id: co2_soilproperties.hh 1327 2009-03-02 14:26:47Z lauser $
#ifndef CO2_SOIL
#define CO2_SOIL

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class Soil: public HomogeneousSoil<Grid, Scalar>
{
public:
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype CoordScalar;
    typedef BasicDomain<Grid, Scalar> ParentType;
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits DomainTraits;
    typedef typename DomainTraits::LocalPosition LocalPosition;
    typedef typename DomainTraits::GlobalPosition GlobalPosition;

    enum
        {    dim=Grid::dimension, numEq=2};

    const FieldMatrix<CoordScalar,dim,dim>& K(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos)
    {
        return permloc_;
    }

    double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.15;
    }

    double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.;
    }

    double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.;
    }

    virtual double heatCap(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return     (800 /* spec. heat cap. of sediment */
                    * 2650 /* density of sediment */
                    * (1-porosity(globalPos, element, localPos)));
    }

    virtual double heatCond(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double sat) const
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

    std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        std::vector<double> param(2);

        //Brooks-Corey parameters
        param[0] = 2.; // lambda
        param[1] = 10000.; // entry-pressure

        return param;
    }

    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid,Scalar>::brooks_corey;
    }


    Soil()
        :HomogeneousSoil<Grid,Scalar>()
    {
        permloc_ = 0;
        permlocWell_ = 0;
        for (int i = 0; i < dim; i++)
            permloc_[i][i] = 2.0e-14;
        for (int i = 0; i < dim; i++)
            permlocWell_[i][i] = 1.0e-12;
    }

private:
    Dune::FieldMatrix<CoordScalar,dim,dim> permloc_, permlocWell_;
};
} // end namespace
#endif
