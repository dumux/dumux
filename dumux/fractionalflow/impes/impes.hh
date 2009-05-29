// $Id$

#ifndef DUNE_IMPES_HH
#define DUNE_IMPES_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, Markus Wolff
 */

namespace Dune
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * coupled diffusion/transport problems
 */

template<class GridView, class Diffusion, class Transport, class VC> class IMPES: public FractionalFlow<
        GridView, Diffusion, Transport, VC>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
typedef    typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FractionalFlow<GridView, Diffusion, Transport, VC> FractionalFlow;
    typedef typename FractionalFlow::RepresentationType PressType;
    typedef typename FractionalFlow::Scalar Scalar;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef typename Transport::RepresentationType RepresentationType;

    virtual void initial()
    {
        Scalar t = 0;
        //initial saturations
        this->initialTransport();
        //call function with true to get a first initialisation of the pressure field
        this->pressure(true,t);
        this->calculateVelocity(t);

        return;
    }

    //calculate saturation defect
    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec,
            Scalar cFLFactor = 1)
    {
        PressType pressOldIter(this->transProblem.variables().pressure());
        PressType pressHelp(this->transProblem.variables().gridSizeDiffusion());
        int satSize = this->transProblem.variables().gridSizeTransport();
        RepresentationType saturation(this->transProblem.variables().saturation());
        RepresentationType satOldIter(this->transProblem.variables().saturation());
        RepresentationType satHelp(satSize);
        RepresentationType satDiff(satSize);
        RepresentationType updateOldIter(satSize);
        RepresentationType updateHelp(satSize);
        RepresentationType updateDiff(satSize);

        //update constitutive functions
        updateMaterialLaws();

        bool converg = false;
        int iter = 0;
        int iterTot = 0;
        updateOldIter = 0;
        while (!converg)
        {
            iter++;
            iterTot++;

            // update pressure: give false as the pressure field is already initialised
            this->pressure(false , t);
            this->calculateVelocity(t);

            //calculate saturation defect
            Transport::update(t, dt, updateVec,cFLFactor, true);

            if (iterFlag)
            { // only needed if iteration has to be done
                this->transProblem.variables().pressure() *= omega;
                pressHelp = pressOldIter;
                pressHelp *= (1-omega);
                this->transProblem.variables().pressure() += pressHelp;
                pressOldIter = this->transProblem.variables().pressure();

                updateHelp = updateVec;
                saturation = this->transProblem.variables().saturation();
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
        this->transProblem.variables().vtkout(name, k);
        return;
    }

    //constitutive functions are updated once if new saturations are calculated and stored in the this->transProblem.variables() object
    void updateMaterialLaws()
    {
        // iterate through leaf grid an evaluate c0 at cell center
        ElementIterator eItEnd = Transport::gridView.template end<0>();
        for (ElementIterator eIt = Transport::gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get cell center in reference element
            const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

            // get global coordinate of cell center
            GlobalPosition globalPos = eIt->geometry().global(localPos);

            int globalIdx = this->transProblem.variables().indexTransport(*eIt);

            Scalar sat = this->transProblem.variables().saturation()[globalIdx];

            std::vector<Scalar> mobilities = this->transProblem.materialLaw().mob(sat, globalPos, *eIt, localPos);

            // initialize mobilities
            this->transProblem.variables().mobilityWetting()[globalIdx]= mobilities[0];
            this->transProblem.variables().mobilityNonWetting()[globalIdx]= mobilities[1];
            this->transProblem.variables().capillaryPressure()[globalIdx]= this->transProblem.materialLaw().pC(sat, globalPos, *eIt, localPos);
            this->transProblem.variables().fracFlowFuncWetting()[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
            this->transProblem.variables().fracFlowFuncNonWetting()[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
        }
        return;
    }

    //! Construct an IMPES object.
    IMPES(Diffusion& diffusion, Transport& transport, int flag = 0, int nIt = 2,
            Scalar maxDef = 1e-5, Scalar om = 1) :
    FractionalFlow(diffusion, transport),
    iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om)
    {
    }

protected:
    const int iterFlag;
    const int nIter;
    const Scalar maxDefect;
    const Scalar omega;
};
}
#endif
