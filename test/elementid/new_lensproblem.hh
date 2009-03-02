/*****************************************************************************
*   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
*   Institute of Hydraulic Engineering                                      *
*   University of Stuttgart, Germany                                        *
*   email: and _at_ poware.org                                              *
*                                                                           *
*   This program is free software; you can redistribute it and/or modify    *
*   it under the terms of the GNU General Public License as published by    *
*   the Free Software Foundation; either version 2 of the License, or       *
*   (at your option) any later version, as long as this copyright notice    *
*   is included in its original form.                                       *
*                                                                           *
*   This program is distributed WITHOUT ANY WARRANTY.                       *
*****************************************************************************/
/*!
* \file LensProblem.hh A Cube of fine sand embedded into a cube of
*                      coarse sand.
*/
#ifndef DUMUX_LENSPROBLEM_HH
#define DUMUX_LENSPROBLEM_HH

#include "lensdomain.hh"
#include "lensnewtoncontroller.hh"

#include <dumux/new_models/2p/pwsnboxmodel.hh>
#include <dumux/new_material/parkerlenhard.hh>
#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/timedisc/new_impliciteulerstep.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <iostream>

namespace Dune
{
namespace Lens
{
// The problem controller class

/** \todo Please doc me! */

template<class ScalarT = double>
class PwSnLensProblem : public PwSnLensDomain<ScalarT>
{
    typedef PwSnLensDomain<ScalarT>         ParentType;
    typedef PwSnLensProblem<ScalarT>        ThisType;
    typedef PwSnBoxModel<ThisType>          Model;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the Pw-Sn model
    typedef typename Model::PwSnTraits          PwSnTraits;
    // the traits for the material relations
    typedef typename ParentType::MaterialTraits MaterialTraits;

private:
    // some constants from the traits for convenience
    enum {
        numEq = BoxTraits::numEq,
        pWIdx = PwSnTraits::pWIdx,
        snIdx = PwSnTraits::snIdx
    };

    // copy some types from the traits for convenience
    typedef typename DomainTraits::Scalar                     Scalar;
    typedef typename DomainTraits::Grid                       Grid;
    typedef typename DomainTraits::Element                    Element;
    typedef typename DomainTraits::ElementIterator            ElementIterator;
    typedef typename DomainTraits::ReferenceElement           ReferenceElement;
    typedef typename DomainTraits::Vertex                     Vertex;
    typedef typename DomainTraits::VertexIterator             VertexIterator;
    typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
    typedef typename DomainTraits::LocalPosition              LocalPosition;
    typedef typename DomainTraits::GlobalPosition             GlobalPosition;

    typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction               SpatialFunction;
    typedef typename BoxTraits::SolutionVector                SolutionVector;
    typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

    typedef typename MaterialTraits::ParkerLenhard            ParkerLenhard;

    typedef Dune::TimeManager                  TimeManager;
    typedef Dune::NewImplicitEulerStep<ThisType>        TimeIntegration;

    typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

    typedef typename Model::NewtonMethod                    NewtonMethod;
    typedef LensNewtonController<NewtonMethod, ThisType>    NewtonController;

public:
    PwSnLensProblem(Scalar initialTimeStepSize, Scalar endTime, const char *dgfFileName)
    : timeManager_(this->grid().comm().rank() == 0),
    model_(*this),
    newtonMethod_(model_),
    newtonCtl_(*this),
    resultWriter_("lens")
    {
        Api::require<Api::BasicDomainTraits, DomainTraits>();
        initGrid_();

        initialTimeStepSize_ = initialTimeStepSize;
        endTime_ = endTime;
    };

    ~PwSnLensProblem()
    {
    }

    //! Actually run the simulation. We outsource the actual work
    //! to the TimeManager here.
    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    }

    ///////////////////////////////////
    // Strings pulled by the TimeManager during the course of the
    // simulation
    ///////////////////////////////////

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // start with a drainage for 30 ksec
        timeManager_.startNextEpisode(DrainEpisode, 30e3);
        timeManager_.setStepSize(initialTimeStepSize_);

        // set the initial condition
        model_.initial();

        // write the inital solution to disk
        writeCurrentResult_();
    }


    /*!
    * \brief Called by the TimeManager in order to get a time
    *        integration on the model.
    *
    * \note timeStepSize and nextStepSize are references and may
    *       be modified by the TimeIntegration. On exit of this
    *       function 'timeStepSize' must contain the step size
    *       actually used by the time integration for the current
    *       steo, and 'nextStepSize' must contain the suggested
    *       step size for the next time step.
    */
    void timeIntegration(Scalar &stepSize, Scalar &nextStepSize)
    {
        // execute the time integration (i.e. Runge-Kutta
        // or Euler).  TODO/FIXME: Note that the time
        // integration modifies the curTimeStepSize_ if it
        // thinks that's appropriate. (IMHO, this is an
        // incorrect abstraction, the simulation
        // controller is responsible for adapting the time
        // step size!)
        timeIntegration_.execute(*this,
                timeManager_.time(),
                stepSize,
                nextStepSize,
                1e100, // firstDt or maxDt, TODO: WTF?
                        endTime_,
                        1.0); // CFL factor (not relevant since we use implicit euler)
    };


    //! called by the TimeIntegration::execute function to let
    //! time pass.
    void updateModel(Scalar &dt, Scalar &nextDt)
    {
        model_.update(dt, nextDt, newtonMethod_, newtonCtl_);
    }

    //! called by the TimeManager whenever a solution for a
    //! timestep has been computed
    void timestepDone()
    {
        // write the current result to disk
        writeCurrentResult_();

        // update the domain with the current solution
        updateDomain_();

        // change the episode of the simulation if necessary
        updateEpisode_();
    };
    ///////////////////////////////////
    // End of simulation control stuff
    ///////////////////////////////////


    ///////////////////////////////////
    // Strings pulled by the PwSnBoxModel during the course of the
    // simulation (-> boundary conditions, initial conditions,
    // etc)
    ///////////////////////////////////

    //! Returns the current time step size in seconds
    Scalar timeStepSize() const
    { return timeManager_.stepSize(); }

    //! Set the time step size in seconds.
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setStepSize(dt); }

    //! evaluate the initial condition for a vert
    void initial(SolutionVector          &values,
            const Element           &element,
            const FVElementGeometry &fvElemGeom,
            int                      scvIdx) const
            {
        const GlobalPosition &globalPos
        = fvElemGeom.subContVol[scvIdx].global;
        /*                const LocalPosition &localPos
                    = fvElemGeom.subContVol[scvIdx].local;
         */

        Scalar pw = -ParentType::densityW() *
        ParentType::gravity()[1] *
        (ParentType::height() - globalPos[1]);

        if (ParentType::onLeftBoundary(globalPos)) {
            Scalar a = -(1 + 0.5/ParentType::height());
            Scalar b = -a*ParentType::upperRight()[1];
            pw = -ParentType::densityW()*
            ParentType::gravity()[1]*(a*globalPos[1] + b);
        }

        Scalar Sn;
        Sn = 0;

        values[pWIdx] = pw; // pw
        values[snIdx] = Sn; // Sn
            }


    // Returns the type of an boundary condition for the wetting
    // phase pressure at a element face
    void boundaryTypes(BoundaryTypeVector         &values,
            const Element              &element,
            const FVElementGeometry    &fvElemGeom,
            const IntersectionIterator &isIt,
            int                         scvIdx,
            int                         boundaryFaceIdx) const
            {
        const GlobalPosition &globalPos
        = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
        //                const LocalPosition &localPos
        //                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

        if (ParentType::onLeftBoundary(globalPos) ||
                ParentType::onRightBoundary(globalPos))
        {
            values[0] = values[1] = Dune::BoundaryConditions::dirichlet;



        }
        else {
            // upper or lower boundary of the grid

            values[0] = values[1] = Dune::BoundaryConditions::neumann;
        }
            }

    //! Evaluate a neumann boundary condition
    void neumann(SolutionVector             &values,
            const Element              &element,
            const FVElementGeometry    &fvElemGeom,
            const IntersectionIterator &isIt,
            int                         scvIdx,
            int                         boundaryFaceIdx) const
            {
        const GlobalPosition &globalPos
        = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
        //                const LocalPosition &localPos
        //                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

        values[pWIdx] = 0;
        values[snIdx] = 0;

        if (ParentType::onUpperBoundary(globalPos)) {
            Scalar relPosX = (ParentType::upperRight()[0] - globalPos[0])/ParentType::width();
            if (0.5 < relPosX && relPosX < 2.0/3.0)
            {
                values[snIdx] = -0.04;
            }
        }

            }


    //! Evaluate a dirichlet boundary condition at a vertex within
    //! an element's face
    void dirichlet(SolutionVector             &values,
            const Element              &element,
            const FVElementGeometry    &fvElemGeom,
            const IntersectionIterator &isIt,
            int                         scvIdx,
            int                         boundaryFaceIdx) const
            {
        const GlobalPosition &globalPos
        = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
        //                const LocalPosition &localPos
        //                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

        Scalar a, b;

        if (ParentType::onLeftBoundary(globalPos))
        {
            a = -(1 + 0.5/ParentType::height());
            b = -a*ParentType::upperRight()[1];
        }
        else {
            a = -1;
            b = ParentType::upperRight()[1];
        }

        values[pWIdx] = -ParentType::densityW()*
        ParentType::gravity()[1]*
        (a*globalPos[1] + b);
        values[snIdx] = 0;
            }

    //! evaluate the mass injection rate of the fluids for a BOX
    //! sub control volume
    void source(SolutionVector          &values,
            const Element           &element,
            const FVElementGeometry &dualElement,
            int                      scvIdx)
    {
        values[pWIdx] = values[snIdx] = 0;
    }

    ///////////////////////////////////
    // End of problem specific stuff
    ///////////////////////////////////

    const FieldVector &gravity() const
         {
            return gravity_;
         }

    //! Return the capillary pressure for a given vert of a element
    Scalar pC(const Element &element,
              int elementIdx,
              int localVertIdx,
              int globalVertIdx,
              Scalar Sw) const
        {

        }

    // return the capillary pressure for a given element
    Scalar porosity(const Element &element) const
        {
        }

    // return the density of the wetting phase
    Scalar densityW() const
        {
        }

    // return the density of the non-wetting phase
    Scalar densityN() const
        {
        }

    // return the density of a phase given by an index. (lower
    // indices means that the phase is more wetting)
    Scalar density(int phase) const
        {
            return (phase == 0)? densityW() : densityN();
        }

    // return the viscosity of the wetting phase
    Scalar viscosityW() const
        {
        }

    // return the viscosity of the non-wetting phase
    Scalar viscosityN() const
        {
        }

    // return the viscosity of a phase given by an index. (lower
    // indices means that the phase is more wetting)
    Scalar viscosity(int phase) const
        {
            return (phase == 0)? viscosityW() : viscosityN();
        }

    // return the mobility of the wetting phase at a vert
    Scalar mobilityW(const Element &element,
                     int elementIdx,
                     int localVertIdx,
                     int globalVertIdx,
                     Scalar Sw) const
        {
        }

    // return the mobility of the non-wetting phase at a vert
    Scalar mobilityN(const Element &element,
                     int elementIdx,
                     int localVertIdx,
                     int globalVertIdx,
                     Scalar Sn) const
        {
        }


private:
    // write results to the output files
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                ParentType::grid().leafView());
        writeVertexFields_(resultWriter_, model_.currentSolution());
        writeElementFields_(resultWriter_, model_.currentSolution());
        resultWriter_.endTimestep();
    }

    // write the fields current solution into an VTK output
    // file.
    void writeVertexFields_(VtkMultiWriter &writer, SpatialFunction &u)
    {
        writer.addScalarVertexFunction("Sn",
                u,
                snIdx);
        writer.addScalarVertexFunction("Pw",
                u,
                pWIdx);

    }


    // simulated time control stuff
    TimeManager     timeManager_;
    TimeIntegration timeIntegration_;
    Scalar          initialTimeStepSize_;
    Scalar          endTime_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;
    VtkMultiWriter   resultWriter_;
};
}
}


#endif
