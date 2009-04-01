// $Id$

#ifndef MATRIXPROPERTIES
#define MATRIXPROPERTIES

#include <dumux/material/property_baseclasses.hh>
#include "dumux/material/randompermeability.hh"

/**
 * \ingroup properties
 * \author Jochen Fritz
 */

namespace Dune
{

/** \todo Please doc me! */

template<class G, class RT>
class HomogeneousSoil: public Matrix2p<G,RT>
{
public:
    typedef    typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1};

    virtual const FieldMatrix<DT,n,n> &K(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi)
    {
        return K_;
    }
    virtual double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 0.3;
    }

    virtual double Sr_w(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
    {
        return 0.;
    }

    virtual double Sr_n(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
    {
        return 0.;
    }

    /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
     * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
    virtual double heatCap(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 790 /* spec. heat cap. of granite */
            * 2700 /* density of granite */
            * (1 - porosity(x, e, xi));
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

    virtual std::vector<double> paramRelPerm(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T) const
    {
        // example for Brooks-Corey parameters
        std::vector<double> param(2);
        param[0] = 2.; // lambda
        param[1] = 0.; // entry-pressures

        //        if (x[0] > 150)
        //            param[0] = 0.5;

        return param;
    }

    virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return Matrix2p<G,RT>::brooks_corey;
    }

    HomogeneousSoil(const double k = 1e-12):Matrix2p<G,RT>()
    {
        K_ = 0;
        for(int i = 0; i < n; i++)
            K_[i][i] = k;
    }

    HomogeneousSoil(const double k1, const double k2, const double k3 = 0, const double k4 = 0, const double k5 = 0, const double k6 = 0, const double k7 = 0, const double k8 = 0, const double k9 = 0):Matrix2p<G,RT>()
    {
        K_ = 0;

        if (n == 1)
        {
            K_[0][0] = k1;
        }
        else if (n == 2)
        {
            K_[0][0] = k1;
            K_[1][1] = k2;
            K_[0][1] = k3;
            K_[1][0] = k4;
        }
        else if (n == 3)
        {
            K_[0][0] = k1;
            K_[1][1] = k2;
            K_[2][2] = k3;
            K_[0][1] = k4;
            K_[0][2] = k5;
            K_[1][0] = k6;
            K_[1][2] = k7;
            K_[2][0] = k8;
            K_[2][1] = k9;     
        }
        else
        {
            std::cout << "Dimension " << n << " is not considered." << std::endl;
        }
    }

    ~HomogeneousSoil()
    {}

private:
    FieldMatrix<DT,n,n> K_;

};

/** \todo Please doc me! */

template<class G, class RT>
class HeterogeneousSoil: public Matrix2p<G,RT>
{
public:
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1};

    const virtual FieldMatrix<DT,n,n> &K (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi)
    {
        return permeability.K(e);
    }

    virtual double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 0.3;
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
            * (1-porosity(x, e, xi));
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
        // example for Brooks-Corey parameters
        std::vector<double> param(2);
        param[0] = 2.; // lambda
        param[1] = 0.; // entry-pressure

        return param;
    }

    virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return Matrix2p<G,RT>::brooks_corey;
    }

    virtual bool readPropertiesFlag()
    {
        return false;
    }

    HeterogeneousSoil(const G& g,const char* name = "permeab.dat", const bool create = true)
        :Matrix2p<G,RT>(),permeability(g, name, create)
    {}

    ~HeterogeneousSoil()
    {}

private:
public:
    RandomPermeability<G> permeability;

};


// Soil with heterogeneous permeability field gained from random simultaion by gstat.
template<class G, class RT>
class GstatHeterogeneousSoil: public Matrix2p<G,RT>
{
public:
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1};

    virtual FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi)
    {
          return permeability.K(e);
    }
    virtual double porosity(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 0.3;
    }

    virtual double Sr_w(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        return 0.2;
    }

    virtual double Sr_n(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        return 0.15;
    }

    /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
             * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
    virtual double heatCap(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return  790 /* spec. heat cap. of granite */
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
        // example for Brooks-Corey parameters
        std::vector<double> param(2);
        param[0] = 2.; // lambda
        param[1] = 0.; // entry-pressure

        return param;
    }

    virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return Matrix2p<G,RT>::linear;
    }

    GstatHeterogeneousSoil(const G& g, const char* gstatCon = "gstatControl.txt", const char* gstatIn = "gstatInput.txt", const char* gstatOut = "permeab.dat", const bool create = true)
        :Matrix2p<G,RT>(),permeability(g, create, gstatOut, gstatCon, gstatIn)
    {}

    ~GstatHeterogeneousSoil()
    {}

private:
public:
    GstatRandomPermeability<G> permeability;
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
