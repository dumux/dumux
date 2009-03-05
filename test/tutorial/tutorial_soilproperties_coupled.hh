#ifndef TUTORIAL_SOILPROPERTIES
#define TUTORIAL_SOILPROPERTIES

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class G, class RT>
class TutorialSoil: public HomogeneousSoil<G, RT> /*@\label{tutorial-coupled:tutorialsoil}@*/
{
public:
    typedef    typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype DT;

    enum{n=G::dimension};

    // function returning the intrinsic permeability tensor K
    // depending on the position within the domain
    const FieldMatrix<DT,n,n> &K (const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-coupled:permeability}@*/
                                  const FieldVector<DT,n>& xi)
    {
        return K_;
    }

    // function returning the porosity of the porous matrix
    // depending on the position within the domain
    double porosity(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-coupled:porosity}@*/
                    const FieldVector<DT,n>& xi) const
    {
        return 0.2;
    }

    // function returning the residual saturation of the wetting fluid
    // depending on the position within the domain and on the temperature
    double Sr_w(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-coupled:srw}@*/
                const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        return 0;
    }

    // function returning the residual saturation of the non-wetting fluid
    // depending on the position within the domain and on the temperature
    double Sr_n(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-coupled:srn}@*/
                const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        return 0;
    }

    // function returning the parameters of the capillary pressure
    // and the relative permeability functions
    // depending on the position within the domain and on the temperature
    std::vector<double> paramRelPerm(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-coupled:parameters}@*/
                                     const FieldVector<DT,n>& xi, const double T = 283.15) const
    {
        std::vector<double> param(2);

        //linear law parameters
        param[0] = 0; // minimal capillary pressure
        param[1] = 0; // maximal capillary pressure

        //Brooks-Corey parameters
        //        param[0] = 2; // lambda
        //        param[1] = 0.; // entry-pressure

        return param;
    }

    // function returning the kind of relation used for the calculation of the capillary
    // pressure and the relative permeabilities depending on the position within the domain
    typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-coupled:flags}@*/
                                                   const FieldVector<DT,n>& xi) const
    {
        return Matrix2p<G,RT>::linear; //flag types defined in
    }                                   //dumux/material/property_baseclasses.hh

    TutorialSoil()
        :HomogeneousSoil<G,RT>(), K_(0)
    {
        for(int i = 0; i < n; i++)
            K_[i][i] = 1e-7;
    }

private:
    FieldMatrix<DT,n,n> K_;
};
} // end namespace
#endif
