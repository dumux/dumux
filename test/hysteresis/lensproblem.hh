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

#include <dumux/new_models/pwsn/pwsnboxmodel.hh>
#include <dumux/new_material/parkerlenhard.hh>
#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/timedisc/new_impliciteulerstep.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <iostream>

#define LENS_WRITE_NEWTON_STEPS 0

namespace Dune
{
namespace Lens
{
    // The problem controller class
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

        // episode control stuff
        enum Episode {
            ImbibEpisode,  // an episode where imbibition of water takes place
            DrainEpisode,  // an episode where drainage of water takes place
            WaitEpisode    // an episode with neither drainage nor imbibition
        };

        typedef Dune::TimeManager<Episode>                  TimeManager;
        typedef Dune::NewImplicitEulerStep<ThisType>        TimeIntegration;

        typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

        typedef typename Model::NewtonMethod                    NewtonMethod;
        typedef LensNewtonController<NewtonMethod, ThisType>    NewtonController;

    public:
        PwSnLensProblem(Scalar initialTimeStepSize, Scalar endTime)
            : timeManager_(this->grid().comm().rank() == 0),
              model_(*this),
              newtonMethod_(model_),
              newtonCtl_(*this),
              resultWriter_("lens")
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

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
        void initial(SolutionVector &dest,
                     const Element &element,
                     GlobalPosition pos,
                     LocalPosition posLocal)
            {
/*                GlobalPosition pos;
                ParentType::vertPosition(pos,
                                            ParentType::vert(element, localVertIdx));
*/
                Scalar pw = -ParentType::densityW() *
                              ParentType::gravity()[1] *
                              (ParentType::height() - pos[1]);

                if (ParentType::onLeftBoundary(pos)) {
                    Scalar a = -(1 + 0.5/ParentType::height());
                    Scalar b = -a*ParentType::upperRight()[1];
                    pw = -ParentType::densityW()*
                          ParentType::gravity()[1]*(a*pos[1] + b);
                }

                Scalar Sn;
#if 0
                if (ParentType::isInLens(pos))
                    Sn = lensMedium_->Snr();
                else
                    Sn = outerMedium_->Snr();
#else
                Sn = 0;
#endif

                dest[0] = pw; // pw
                dest[1] = Sn; // Sn
            }


        // Returns the type of an boundary contition for the wetting
        // phase pressure at a element face
        void boundaryTypes(BoundaryTypeVector &dest,
                           const Element &element,
                           const IntersectionIterator &face,
                           const GlobalPosition &pos,
                           const LocalPosition &localPos)

            {
                // get the integration point of the boundary face in
                // world coodinates
//        GlobalPosition &pos = dualElement.boundaryFace[dcBFIdx].ipGlobal;

                if (ParentType::onLeftBoundary(pos) ||
                    ParentType::onRightBoundary(pos))
                {
                    dest[0] = dest[1] = Dune::BoundaryConditions::dirichlet;

#if !USE_ORIG_PROB
                    Scalar relPosY = (pos[1] - ParentType::lowerLeft()[1])/ParentType::height();
                    if (relPosY < 0.80)
                    {
                        dest[0] = dest[1] = Dune::BoundaryConditions::neumann;
                    }
#endif


                }
                else {
                    // upper or lower boundary of the grid

#if 0
                    if (ParentType::onUpperBoundary(pos))
                        dest[0] = dest[1] = Dune::BoundaryConditions::dirichlet;
                    else
                        dest[0] = dest[1] = Dune::BoundaryConditions::neumann;
#else
                    dest[0] = dest[1] = Dune::BoundaryConditions::neumann;
#endif
                }
            }

        //! Evaluate a neumann boundary condition
        void neumann(SolutionVector &dest,
                     const Element &element,
                     const IntersectionIterator &face,
                     const GlobalPosition &pos,
                     const LocalPosition &localPos)
            {
                dest[pWIdx] = 0;
                dest[snIdx] = 0;

                // get the integration point of the boundary face in
                // world coodinates
//        GlobalPosition &pos = dualElement.boundaryFace[dcBFIdx].ipGlobal;

                if (ParentType::onUpperBoundary(pos)) {
#if USE_ORIG_PROB
                    Scalar relPosX = (ParentType::upperRight()[0] - pos[0])/ParentType::width();
                    if (0.5 < relPosX && relPosX < 2.0/3.0)
                    {
                        dest[snIdx] = -0.04;
                    }
#endif
                }
#if !USE_ORIG_PROB
                else if (ParentType::onLowerBoundary(pos)) {
                    dest[snIdx] = 0.0;
                    if (timeManager_.episode() == DrainEpisode)
                        // drain water
                        dest[snIdx] = -0.04;
                    else if (timeManager_.episode() == ImbibEpisode)
                        // imbibition of water
                        dest[snIdx] = 0.04;
                }
#endif

            }


        //! Evaluate a dirichlet boundary condition at a vert within
        //! an element's face
        void dirichlet(SolutionVector &dest,
                       const Element &element,
                       int   vertIdx,
                       int   globalVertexIdx)

/*        void dirichlet(SolutionVector &dest,
                       const Element &element,
                       const IntersectionIterator &face,
                       const GlobalPosition &pos,
                       const LocalPosition &localPos)
*/
            {
                Scalar a, b;

                const LocalPosition &localPos = element.geometry().corner(vertIdx);
                GlobalPosition pos = element.geometry().global(localPos);

                if (ParentType::onLeftBoundary(pos))
                {
                    a = -(1 + 0.5/ParentType::height());
                    b = -a*ParentType::upperRight()[1];
                }
                else {
                    a = -1;
                    b = ParentType::upperRight()[1];
                }

                dest[pWIdx] = -ParentType::densityW()*
                                  ParentType::gravity()[1]*
                                  (a*pos[1] + b);
#if 0
                dest[snIdx] = ParentType::outerMedium().Snr();
#else
                dest[snIdx] = 0;
#endif
            }

        //! evaluate the mass injection rate of the fluids for a BOX
        //! sub control volume
        void source(SolutionVector &dest,
                        const Element &element,
                        const FVElementGeometry &dualElement,
                        int subControlVolumeId)
            {
                dest[pWIdx] = dest[snIdx] = 0;
            }

        ///////////////////////////////////
        // End of problem specific stuff
        ///////////////////////////////////

        //! called by the LensNewtonContoller when the newton method
        //! is started.
        void newtonBegin()
            {
#if LENS_WRITE_NEWTON_STEPS
                convergenceWriter_ =
                    new VtkMultiWriter((boost::format("lens-convergence-t=%.2f-dt=%.2f")
                                        %timeManager_.time()
                                        %model_.localJacobian().getDt()).str());
#endif // LENS_WRITE_NEWTON_STEPS
            }

        //! called by the LensNewtonContoller when a newton step is
        //! finished.
        void newtonEndStep(SpatialFunction &u, SpatialFunction &uOld)
            {
#if LENS_WRITE_NEWTON_STEPS
                if (newtonCtl_.newtonNumSteps() == 1) {
                    convergenceWriter_->beginTimestep(0,
                                                      ParentType::grid().leafView());
                    writeConvergenceFields_(uOld, uOld);
                    convergenceWriter_->endTimestep();
                }


                convergenceWriter_->beginTimestep(newtonCtl_.newtonNumSteps(),
                                                  ParentType::grid().leafView());
                writeConvergenceFields_(u, uOld);
                convergenceWriter_->endTimestep();
#endif // LENS_WRITE_NEWTON_STEPS
            }

        //! called by the LensNewtonContoller when the newton method
        //! is finished.
        void newtonEnd()
            {
#if LENS_WRITE_NEWTON_STEPS
                delete convergenceWriter_;
#endif // LENS_WRITE_NEWTON_STEPS
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

#if LENS_WRITE_NEWTON_STEPS
        void writeConvergenceFields_(SpatialFunction &u, SpatialFunction &uOld)
            {
                SpatialFunction diff(ParentType::grid());
                for (int i=0; i < (*diff).size(); ++i) {
                    (*diff)[i] = (*u)[i] - (*uOld)[i];
                    if ((*diff)[i][0] > 1e12)
                        (*diff)[i][0] = 1e12;
                    if ((*diff)[i][1] > 1e12)
                        (*diff)[i][1] = 1e12;
                }

                convergenceWriter_->addScalarVertexFunction("reduction Sn",
                                                            diff,
                                                            ParentType::vertMap(),
                                                            snIdx);
                convergenceWriter_->addScalarVertexFunction("reduction Pw",
                                                            diff,
                                                            ParentType::vertMap(),
                                                            pWIdx);
                writeVertexFields_(*convergenceWriter_, u);
                writeElementFields_(*convergenceWriter_, u);
            };
#endif // LENS_WRITE_NEWTON_STEPS


        // appends a element centered capillary pressure field to the
        // multi writer's current timestep
        void writeElementFields_(VtkMultiWriter &writer, SpatialFunction &u)
            {
                // TODO
                /*
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
                int nElements =  ParentType::numElements();
                ScalarField *pC = writer.template createField<Scalar, 1>(nElements);
                ScalarField *dpC_dSw = writer.template createField<Scalar, 1>(nElements);

                ElementIterator it = ParentType::grid().template leafbegin<0>();
                ElementIterator endIt = ParentType::grid().template leafend<0>();
                for (; it != endIt; ++it) {
                    // extract the current solution's Sn component
                    const ElementReferenceElement &refElem =
                        DomainTraits::referenceElement(it->geometry().type());
                    Scalar Sn = u.evallocal(snIdx,
                                             *it,
                                             refElem.position(0,0));

                    // calculate the capillary pressure
                    ElementState &eState = ParentType::elementState(*it);
                    int eIdx = ParentType::elementIdx(*it);
                    Scalar Sw = 1 - Sn;
                    (*pC)[eIdx] = ParentType::pC(eState, Sw);
                    (*dpC_dSw)[eIdx] = ParentType::dpC_dSw(eState, Sw);
                }
                writer.addElementData(pC, "capillary pressure");
                writer.addElementData(dpC_dSw, "dpC/dSw");
                */
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

                /*
                SpatialFunction globResidual(ParentType::grid());
                model_.evalGlobalResidual(globResidual);
                writer.addScalarVertexFunction("global residual Sn",
                                               globResidual,
                                               snIdx);
                writer.addScalarVertexFunction("global residual Pw",
                                               globResidual,
                                               pWIdx);
                */
            }

        // called whenever a solution for a timestep has been computed.
        void updateDomain_()
            {
#if USE_HYSTERESIS
#if USE_NODE_PARAMETERS
                VertexIterator it = ParentType::vertexBegin();
                VertexIterator endit = ParentType::vertexEnd();
                for (; it != endit; ++it) {
                    int vertIdx = ParentType::vertIdx(*it);
                    Scalar Sn = (*model_.currentSolution())[vertIdx][snIdx];
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(vertexState(*it), Sw);
                }
#else
                // update the parker-lenhard hystersis model
                // for all elements
                ElementIterator it = ParentType::grid().template leafbegin<0>();
                ElementIterator endit = ParentType::grid().template leafend<0>();
                for (; it != endit; ++it) {
                    // get the barycenter of the current element
                    const Dune::GeometryType &geoType = it->type();
                    const LocalPosition &localPos =
                        DomainTraits::referenceElement(geoType).position(0,0);

                    // evaluate the solution for the non-wetting
                    // saturation at this point
                    Scalar Sn = model_.u().evallocal(snIdx,
                                                      *it,
                                                      localPos);
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(elementState(elementIdx(*it)), Sw);
                }
#endif
#endif
            }

        void updateEpisode_()
            {
                Scalar    len = timeManager_.episodeLength();
                Episode   epi = timeManager_.episode();
                int        epiIdx = timeManager_.episodeIndex();

                if (timeManager_.time() >= endTime_) {
                    timeManager_.setFinished();
                    return;
                }

#if !USE_ORIG_PROB
                if (!timeManager_.episodeIsOver())
                    return;


                switch (epiIdx) {
                    case 1:
                        timeManager_.startNextEpisode(WaitEpisode, 50e3);
                        return;
                    case 2:
                        timeManager_.startNextEpisode(ImbibEpisode, 10e3);
                        return;
                    case 3:
                        timeManager_.startNextEpisode(WaitEpisode, 50e3);
                        return;
                    case 4:
                        timeManager_.startNextEpisode(DrainEpisode, 20e3);
                        return;
                }

                switch (epi) {
                    case ImbibEpisode:
                        len *= 2;
                        len = std::max(5e3, len);
                        epi = DrainEpisode;
                        break;
                    case DrainEpisode:
                        len /= 3;
                        len = std::max(2e3, len);
                        epi = ImbibEpisode;
                        break;
                    default:
                        throw "ooops: unexpected episode type";
                }

                timeManager_.startNextEpisode(epi, len);
#endif
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
#if LENS_WRITE_NEWTON_STEPS
        VtkMultiWriter *convergenceWriter_;
#endif // LENS_WRITE_NEWTON_STEPS

    };
}
}


#endif
