// $Id$

#ifndef DUNE_FVSPLITTEDTRANSPORT_HH
#define DUNE_FVSPLITTEDTRANSPORT_HH

#include "dumux/transport/splittedtransport.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * \defgroup transport DeprecatedTransport
 */

namespace Dune
{
//! \ingroup transport
//! Base class for defining an instance of a numerical transport model.
/*! An interface for defining a numerical transport model for the
 *  solution of equations of the form
 *  \f$S_t - \text{div}\, (f_\text{w}(S) \boldsymbol{v}_\text{total}) = 0\f$,
 * \f$S = g\f$ on \f$\Gamma_1\f$, and \f$S(t = 0) = S_0\f$. Here,
 * \f$S\f$ denotes the wetting phase saturation,
 * \f$\boldsymbol{v}_\text{total}\f$ the total velocity,
 * and \f$f_\text{w}\f$ the wetting phase fractional flow function.

 - Grid      a DUNE grid type
 - RT        type used for return values
 - RepresentationType   type of the vector holding the saturation values
 - VelType   type of the vector holding the velocity values

*/
template<class G, class RT, class VC>
class FVSplittedTransport : public SplittedTransport< G, RT, VC,
                                                      DeprecatedFVTransport<G, RT, VC>, DeprecatedFVTransport<G, RT, VC> >
{
public:
    typedef typename DeprecatedFVTransport<G, RT, VC>::RepresentationType HyperbolicRepresentationType;
    typedef typename DeprecatedFVTransport<G, RT, VC>::RepresentationType ParabolicRepresentationType;
    typedef typename DeprecatedFVTransport<G, RT, VC>::RepresentationType RepresentationType;
    typedef DeprecatedFVTransport<G, RT, VC>  HyperbolicType;
    typedef DeprecatedFVTransport<G, RT, VC>  ParabolicType;


    virtual void transferHyperbolicToParabolic(const HyperbolicRepresentationType& hyperSat,
                                               ParabolicRepresentationType& paraSat)
    {
        paraSat = hyperSat;
        return;
    }

    virtual void transferParabolicToHyperbolic(const ParabolicRepresentationType& paraSat,
                                               HyperbolicRepresentationType& hyperSat)
    {
        hyperSat = paraSat;
        return;
    }

    virtual void transferParabolicToRepresentationType(const ParabolicRepresentationType& paraSat,
                                                       RepresentationType& rSat)
    {
        rSat = paraSat;
        return;
    }

    virtual void transferRepresentationTypeToParabolic(const RepresentationType& rSat,
                                                       ParabolicRepresentationType& paraSat)
    {
        paraSat = rSat;
        return;
    }

    //! generate vtk output
    virtual void vtkout (const char* name, int k) const
    {
        VTKWriter<typename G::LevelGridView> vtkwriter(this->grid.levelView(this->parabolicLevel()));

        char fname[128];
        sprintf(fname,"%s-%05d",name,k);
        vtkwriter.addCellData(this->sat,"saturation");
        vtkwriter.write(fname,Dune::VTKOptions::ascii);
    }

    /*! @brief constructor
     *  @param g a DUNE grid object
     *  @param prob an object of class DeprecatedTransportProblem or derived
     */
    FVSplittedTransport(const G& g, HyperbolicType& hyper, ParabolicType& para)
        : SplittedTransport< G, RT, VC, DeprecatedFVTransport<G, RT, VC>, DeprecatedFVTransport<G, RT, VC> >(g, hyper, para)
    {
        this->sat.resize(this->parabolicPart.transproblem.variables.saturation.size());
    }
};

}
#endif
