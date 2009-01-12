/// $Id$
#ifndef CO2_SOIL
#define CO2_SOIL

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class G, class RT>
class Soil: public HomogeneousSoil<G, RT>
{
public:
typedef    typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype DT;
    enum
    {    dim=G::dimension, m=2};

    const FieldMatrix<DT,dim,dim>& K(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
    {
           return permloc_;
    }

    double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return 0.15;
    }

    double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
    {
        return 0.;
    }

    double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
    {
        return 0.;
    }

    virtual double heatCap(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return     (800 /* spec. heat cap. of sediment */
                        * 2650 /* density of sediment */
                        * (1-porosity(x, e, xi)));
    }

    virtual double heatCond(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double sat) const
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

    std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T = 283.15) const
    {
        std::vector<double> param(2);

        //Brooks-Corey parameters
        param[0] = 2.; // lambda
        param[1] = 10000.; // entry-pressure

        return param;
    }

    typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return Matrix2p<G,RT>::brooks_corey;
    }


    Soil()
    :HomogeneousSoil<G,RT>()
    {
      permloc_ = 0;
      permlocWell_ = 0;
      for (int i = 0; i < dim; i++)
        permloc_[i][i] = 2.0e-14;
          for (int i = 0; i < dim; i++)
               permlocWell_[i][i] = 1.0e-12;
    }

    private:
    Dune::FieldMatrix<DT,dim,dim> permloc_, permlocWell_;
};
} // end namespace
#endif
