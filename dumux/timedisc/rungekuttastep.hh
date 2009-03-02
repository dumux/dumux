// $Id$

#ifndef DUNE_RUNGEKUTTASTEP_HH
#define DUNE_RUNGEKUTTASTEP_HH

namespace Dune {

/** \todo Please doc me! */

template<class G, class Model>
class RungeKuttaStep : public TimeStep<G, Model>
{
public:
    void execute(Model& model, double t, double& dt,
                 double maxDt, double tEnd, double cFLFactor)
    {
        // allocate temporary vectors for the updates
        typedef typename Model::RepresentationType RepresentationType;
        RepresentationType k0(*model);
        RepresentationType k1(*model);
        RepresentationType k2(*model);
        RepresentationType k3(*model);
        RepresentationType k4(*model);
        RepresentationType help(*model);
        double dummy;

        // scale dt with safety factor
        dt = std::min( dt, maxDt );
        dt = std::min( dt, tEnd - t);

        // obtain the first update and the time step size
        model.update(t, dt, k1,cFLFactor);

        // scale dt with safety factor
        dt *= cFLFactor;
        dt = std::min( dt, maxDt );
        dt = std::min( dt, tEnd - t);

        double theta = 0.25;
        switch (stages) {
        case 1: // Euler
            // explicit Euler: Sat <- Sat + dt*N(Sat)
            *model += (k1 *= dt);
            break;
        case 2: // Heun
            // second stage: k2 = N(Sat + dt*k1)
            help = k1;
            *model += (help *= dt);
            model.update(t, dummy, k2,cFLFactor);

            // combine Sat <- Sat + 0.5*dt*(k1 + k2);
            *model = k1;
            *model += k2;
            *model *= 0.5*dt;
            *model += k0;
            break;
        case 3:
            // second stage: k2 = N(Sat + dt*k1)
            help = k1;
            *model += (help *= dt);
            model.update(t, dummy, k2,cFLFactor);

            // third stage: k3 = N(Sat + 0.5*dt*(k1 + k2))
            *model = k0;
            help = k1;
            *model += (help *= (1.0 - theta)*dt);
            help = k2;
            *model += (help *= theta*dt);
            model.update(t, dummy, k3,cFLFactor);

            // combine Sat <- Sat + 0.5*dt*(k1 + k3);
            *model = (k1 *= 1.0 -theta);
            *model += (k3 *= theta);
            *model *= dt;
            *model += k0;
            break;
        case 4:
            // second stage: k2 = N(Sat + 0.5*dt*k1)
            help = k1;
            *model += (help *= 0.5*dt);
            model.update(t, dummy, k2,cFLFactor);

            // third stage: k3 = N(Sat + 0.5*dt*k2)
            *model = k0;
            help = k2;
            *model += (help *= 0.5*dt);
            model.update(t, dummy, k3,cFLFactor);

            // fourth stage: k4 = N(Sat + dt*k3)
            *model = k0;
            help = k3;
            *model += (help *= dt);
            model.update(t, dummy, k4,cFLFactor);

            // combine Sat <- Sat + 1/6*dt*(k1 + 2k2 + 2k3 + k4);
            *model = k1;
            *model += (k2 *= 2.0);
            *model += (k3 *= 2.0);
            *model += k4;
            *model *= dt/6.0;
            *model += k0;
            break;
        default:
            std::cout << "number of stages (" << stages << ") not implemented for RK scheme!" << std::endl;
            break;
        }

        return;
    }

    RungeKuttaStep(int s = 1)
        : stages(s)
    { }

private:
    int stages;
};
}
#endif
