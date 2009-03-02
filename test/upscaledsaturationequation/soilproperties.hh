#ifndef MATRIXPROPERTIES
#define MATRIXPROPERTIES

#include <dumux/material/property_baseclasses.hh>
#include "dumux/material/randompermeability.hh"

namespace Dune
{

template<class G, class RT>
class HeterogeneousSoil: public Matrix2p<G, RT>
{
public:
    enum
        {
            n = G::dimension, m = 1
        };

    typedef    typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype DT;
    typedef Dune::FieldMatrix<double, n, n> MB;
    typedef Dune::Matrix<MB> DispersionType;
    typedef Dune::Matrix<FieldVector<double,1> > DispersionSatType;
    typedef Dune::FieldVector<double, 2*n> VB;
    typedef Dune::Matrix<VB> ConvectionType;

    virtual FieldMatrix<DT,n,n> K (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi)
    {
        //        FieldMatrix<DT,n,n> K(0);
        //
        //        for(int i = 0; i < n; i++)
        //        K[i][i] = 1e-10;
        //
        //        return K;
        return permeability.K(e);
    }

    virtual FieldMatrix<DT,n,n> D (double& sat, int index)
    {
        return 0;// dispersion[]
    }

    virtual DispersionType& getDispersion()
    {
        return dispersion;
    }

    virtual DispersionSatType& getDispersionSat()
    {
        return dispersionSat;
    }

    virtual ConvectionType& getM()
    {
        return convectiveCorrection;
    }

    virtual ConvectionType& getMSat()
    {
        return convectiveCorrectionSat;
    }

    virtual double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 0.2;
    }

    virtual double Sr_w(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        return 0;
    }

    virtual double Sr_n(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        return 0;
    }

    /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
     * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
    virtual double heatCap(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 790 /* spec. heat cap. of granite */
            * 2700 /* density of granite */
            * porosity(x, e, xi);
    }

    virtual double heatCond(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double sat) const
    {
        static const double lWater = 0.6;
        static const double lGranite = 2.8;
        double poro = porosity(x, e, xi);
        double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
        double ldry = pow(lGranite, (1-poro));
        return ldry + sqrt(sat) * (ldry - lsat);
    }

    virtual std::vector<double> paramRelPerm(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        std::vector<double> param(2);

        //linear parameters
        param[0] = 0;
        param[1] = 0.;

        //Brooks-Corey parameters
        //                    param[0] = 2; // lambda
        //                    param[1] = 0.; // entry-pressure

        return param;
    }

    typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        //                    return Matrix2p<G,RT>::brooks_corey;
        return Matrix2p<G,RT>::linear;
    }

    HeterogeneousSoil(const G& g,const char* name = "permeab.dat", const bool create = true)
        :Matrix2p<G,RT>(),permeability(g, g.maxLevel(), name, create),dispersion(),dispersionSat()
    {}

    ~HeterogeneousSoil()
    {}

private:
public:
    LevelRandomPermeability<G> permeability;
    DispersionType dispersion;
    DispersionSatType dispersionSat;
    ConvectionType convectiveCorrection;
    ConvectionType convectiveCorrectionSat;
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
