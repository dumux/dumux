// $Id$

/*!
 * \file
 * \brief A 1p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_1P2C_NEWTON_CONTROLLER_HH
#define DUMUX_1P2C_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/new_newtoncontroller.hh>

namespace Dune {
/*!
 * \brief A 1p2c specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class NewtonMethod>
class OnePTwoCNewtonController
    : public NewtonControllerBase<NewtonMethod, OnePTwoCNewtonController<NewtonMethod> >
{
public:
    typedef OnePTwoCNewtonController<NewtonMethod>        ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType>  ParentType;

    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    OnePTwoCNewtonController(Scalar tolerance = 1e-5,
                               int targetSteps = 9,
                               int maxSteps = 18)
        : ParentType(tolerance, targetSteps, maxSteps)
    {};

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    void newtonEndStep(Function &u, Function &uOld)
    {
        // call the method of the base class
        ParentType::newtonEndStep(u, uOld);
    }

    //! Suggest a new time stepsize based either on the number of newton
    //! iterations required or on the variable switch
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        /*
          #warning "HACK: remove this:"
          return std::min(1e6,
          ParentType::suggestTimeStepSize(oldTimeStep));
          // end hack
          */

        // use function of the newtoncontroller
        return ParentType::suggestTimeStepSize(oldTimeStep);
    }

    //! Returns true if another iteration should be done.
    bool newtonProceed(Function &u)
    {
        return ParentType::newtonProceed(u);

        bool baseProceed = ParentType::newtonProceed(u);

        return baseProceed;
    }

    /** \todo Please doc me! */

protected:
    friend class NewtonControllerBase<NewtonMethod, ThisType>;
    //! called by the base class the get an indication of how physical
    //! an iterative solution is 1 means "completely physical", 0 means
    //! "completely unphysical"
    Scalar physicalness_(Function &u)
    {
        return 1.0;

    }
};
}

#endif
