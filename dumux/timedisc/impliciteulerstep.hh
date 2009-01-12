// $Id$

#ifndef DUNE_IMPLICITEULERSTEP_HH
#define DUNE_IMPLICITEULERSTEP_HH

namespace Dune {
    template<class G, class Model>
    class ImplicitEulerStep : public TimeStep<G, Model>
    {
    public:
        void execute(Model& model, double t, double& dt,
                        double maxDt, double tEnd, double cFLFactor)
        {
            double eps = 1e-8;
                        dt = std::min( dt, maxDt );

            if (tEnd - t <= (1+eps)*dt)
                dt = (1+eps)*(tEnd - t);

            model.update(dt);
                dt = std::min(dt, tEnd - t);

            return;
        }

    };
}
#endif



