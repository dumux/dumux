// $Id: impes.hh 1329 2009-03-02 15:14:45Z lauser $

#ifndef DUNE_IMPES2P_HH
#define DUNE_IMPES2P_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow_phasepressure.hh"

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

template<class Grid, class Diffusion, class Transport, class VC> class IMPES: public FractionalFlow<
        Grid, Diffusion, Transport, VC>
{
    enum
    {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };
typedef    typename Grid::LevelGridView GridView;
typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FractionalFlow<Grid, Diffusion, Transport, VC> FractionalFlow;
    typedef typename FractionalFlow::RepresentationType PressType;
    typedef typename FractionalFlow::Scalar Scalar;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef typename Transport::RepresentationType RepresentationType;

    //dummy to overload abstract function in the baseclass - not used any longer!!!
    virtual void totalVelocity(const Scalar t=0)
    {
        return;
    }

    virtual void initial()
    {
        Scalar t = 0;
        //initial saturations
        this->initialTransport();
        //calculate constitutive functions
        updateMaterialLaws();
        //call function with true to get a first initialisation of the pressure field
        this->pressure(true,t);

        return;
    }

    //calculate saturation defect
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

            //calculate saturation defect
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

    //constitutive functions are updated once if new saturations are calculated and stored in the variables object
    void updateMaterialLaws()
    {
        const GridView& gridView = this->grid_.levelView(this->transProblem.variables.levelTransport);
        // iterate through leaf grid an evaluate c0 at cell center
        ElementIterator eItEnd = gridView.template end<0>();
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get cell center in reference element
            const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

            // get global coordinate of cell center
            GlobalPosition globalPos = eIt->geometry().global(localPos);

            int globalIdx = this->transProblem.variables.indexSetTransport.index(*eIt);

            Scalar sat = this->transProblem.variables.saturation[globalIdx];

            std::vector<Scalar> mobilities = this->transProblem.materialLaw().mob(sat, globalPos, *eIt, localPos);

            // initialize mobilities
            this->transProblem.variables.mobilityWetting[globalIdx]= mobilities[0];
            this->transProblem.variables.mobilityNonWetting[globalIdx]= mobilities[1];
            this->transProblem.variables.capillaryPressure[globalIdx]= this->transProblem.materialLaw().pC(sat, globalPos, *eIt, localPos);
            this->transProblem.variables.fracFlowFuncWetting[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
            this->transProblem.variables.fracFlowFuncNonWetting[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
        }
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
