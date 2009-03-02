// $Id$

#ifndef DUNE_IMPES_HH
#define DUNE_IMPES_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, last changed by Markus Wolff
 */

namespace Dune
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * coupled diffusion/transport problems
 */

template<class Grid, class Diffusion, class Transport, class VC> class IMPES: public FractionalFlow<
    Grid, Diffusion, Transport, VC>
{

    typedef Dune::FractionalFlow<Grid, Diffusion, Transport, VC> FractionalFlow;
    typedef    typename FractionalFlow::RepresentationType PressType;
    typedef typename FractionalFlow::Scalar Scalar;

public:
    typedef typename Transport::RepresentationType RepresentationType;

    virtual void totalVelocity(const Scalar t=0)
    {
        this->calcTotalVelocity(t);
        return;
    }

    virtual void initial()
    {
        Scalar t = 0;
        this->initialTransport();
        this->pressure(t);
        totalVelocity(t);
        return;
    }

    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec,
                       Scalar cFLFactor = 1)
    {
        int pressSize = variables.pressure.size();
        PressType pressOldIter(variables.pressure);
        PressType pressHelp(pressSize);
        int satSize = variables.saturation.size();
        RepresentationType saturation(variables.saturation);
        RepresentationType satOldIter(variables.saturation);
        RepresentationType satHelp(satSize);
        RepresentationType satDiff(satSize);
        RepresentationType updateOldIter(satSize);
        RepresentationType updateHelp(satSize);
        RepresentationType updateDiff(satSize);

        bool converg = false;
        int iter = 0;
        int iterTot = 0;
        updateOldIter = 0;
        while (!converg)
        {
            iter++;
            iterTot++;
            // update pressure
            pressure(t);
            totalVelocity(t);

            Transport::update(t, dt, updateVec,cFLFactor);
            if (iterFlag)
            { // only needed if iteration has to be done
                variables.pressure *= omega;
                pressHelp = pressOldIter;
                pressHelp *= (1-omega);
                variables.pressure += pressHelp;
                pressOldIter = variables.pressure;

                updateHelp = updateVec;
                saturation = variables.saturation;
                saturation += (updateHelp *= (dt*cFLFactor));
                saturation *= omega;
                satHelp = satOldIter;
                satHelp *= (1-omega);
                saturation += satHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                satOldIter = saturation;
                updateOldIter = updateVec;
            }
            // break criteria for iteration loop
            if (iterFlag==2&& dt*updateDiff.two_norm()/(saturation).two_norm() <= maxDefect )
                converg = true;
            else if (iterFlag==2 && (saturation.infinity_norm()> 1 || saturation.two_norm()> 1))
            {
                converg = false;
            }
            else if (iterFlag==2&& iter> nIter )
            {
                std::cout << "Nonlinear loop in IMPES.update exceeded nIter = "
                          << nIter << " iterations."<< std::endl;
                return 1;
            }
            else if (iterFlag==1&& iter> nIter )
                converg = true;
            else if (iterFlag==0)
                converg = true;
        }
        // outputs
        if (iterFlag==2)
            std::cout << "Iteration steps: "<< iterTot << std::endl;
        std::cout.setf(std::ios::scientific, std::ios::floatfield);

        return 0;
    }

    virtual void vtkout(const char* name, int k) const
    {
        variables.vtkout(name, k);
        return;
    }

    //! Construct an IMPES object.
    IMPES(Diffusion& diffusion, Transport& transport, int flag = 0, int nIt = 2,
          Scalar maxDef = 1e-5, Scalar om = 1) :
        FractionalFlow(diffusion, transport),
        iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om), variables(this->transProblem.variables)
    {
    }

protected:
    const int iterFlag;
    const int nIter;
    const Scalar maxDefect;
    const Scalar omega;
    VC& variables;
};
}
#endif
