// $Id$

#ifndef DUNE_TIMESTEP_HH
#define DUNE_TIMESTEP_HH

namespace Dune {

/** \todo Please doc me! */

    template<class G, class Model>
    class TimeStep
    {
    public:
        virtual void execute(Model& model, double t, double& dt,
                                double maxDt, double tEnd, double cFLFactor) = 0;

        //! always define virtual destructor in abstract base class
        virtual ~TimeStep () {}
    };
}
#endif
