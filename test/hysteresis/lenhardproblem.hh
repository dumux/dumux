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
 * \file
 *
 * \brief The 1D-Lenhard benchmark problem for hysteresis.
 *
 * This has been modeled to match as closely as possible to the one
 * described at:
 *
 * Sheta, Hussam: "Simulation von Mehrphasenvorgaengen in poroesen
 *     Medien unter Einbeziehung von Hystereseeffekten", PhD theses,
 *     Braunschweig 1999, pp. 112
 */
#ifndef DUMUX_LENHARDPROBLEM_HH
#define DUMUX_LENHARDPROBLEM_HH

#include "lenharddomain.hh"

#include <dumux/new_models/2p/pwsnboxmodel.hh>
#include <dumux/new_models/2p/pwsnnewtoncontroller.hh>

#include <dumux/new_material/parkerlenhard.hh>
#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/timedisc/new_impliciteulerstep.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include "lenhardnewtoncontroller.hh"

#include <iostream>

namespace Dune
{
namespace Lenhard
{
    // The problem controller class
    template<class ScalarT>
    class PwSnLenhardProblem : public PwSnLenhardDomain<ScalarT>
    {
        typedef PwSnLenhardDomain<ScalarT>  ParentType;
        typedef PwSnLenhardProblem<ScalarT> ThisType;
        typedef PwSnBoxModel<ThisType>      Model;

    public:
        //! the traits of the domain
        typedef typename ParentType::DomainTraits   DomainTraits;
        //! the traits of the BOX scheme
        typedef typename Model::BoxTraits           BoxTraits;
        //! the traits of the Pw-Sn model
        typedef typename Model::PwSnTraits          PwSnTraits;
        //! the traits for the material relations
        typedef typename ParentType::MaterialTraits MaterialTraits;

    private:
        // some constants from the traits for convenience
        enum {
            numEq   = BoxTraits::numEq,
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
        typedef typename MaterialTraits::ElementState             ElementState;
        typedef typename MaterialTraits::VertexState              VertexState;

        // episode control stuff
        enum Episode {
            ImbibEpisode,  // an episode where imbibition of water takes place
            DrainEpisode,  // an episode where drainage of water takes place
            WaitEpisode    // an episode with neither drainage nor imbibition
        };

        typedef Dune::TimeManager<Episode>           TimeManager;
        typedef Dune::NewImplicitEulerStep<ThisType> TimeIntegration;
        typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

//        typedef PwSnNewtonController<Model>   NewtonController;
        typedef typename Model::NewtonMethod                    NewtonMethod;
        typedef LenhardNewtonController<NewtonMethod, ThisType> NewtonController;

    public:
        PwSnLenhardProblem(Scalar initialTimeStepSize, Scalar endTime)
            : timeManager_(this->grid().comm().rank() == 0),
              model_(*this),
              newtonMethod_(model_),
              newtonCtl_(*this),
#if LENHARD_EXPERIMENT == 1
              resultWriter_("LenhardExp1")
#elif LENHARD_EXPERIMENT == 2
              resultWriter_("LenhardExp2")
#endif
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

#ifdef USE_NODE_PARAMETERS
                maxPc_ = ParkerLenhard::pC(ParentType::vertexState(0),
                                           0.0);
#else // !USE_NODE_PARAMETERS
                maxPc_ = ParkerLenhard::pC(ParentType::elementState(0),
                                           0.0);
#endif
                initialTimeStepSize_ = initialTimeStepSize;
                endTime_ = endTime;
            };

        ~PwSnLenhardProblem()
            {
            }

        //! Actually run the simulation. We outsource the actual work
        //! to the TimeManager here.
        bool simulate()
            {
                Dune::Timer timer;
                timer.reset();

                timeManager_.runSimulation(*this);

                std::cout <<
                    boost::format("LenhardSimulation took %.3f seconds\n")
                    %timer.elapsed();
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
                // set the initial water level and wait for 10 minutes
                updateEpisode_();

                 // set the initial solution
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
                
                initial_(values, globalPos);
            }

    private:
        // internal method for the initial condition (reused for the
        // dirichlet conditions!)
        void initial_(SolutionVector       &dest,
                      const GlobalPosition &globalPos) const
        {
                Scalar h = curHydraulicHead_ - globalPos[0];
                Scalar pH = hydrostaticPressure_(h);

                // initially the lower 67cm are filled with water the
                // remaining 5 with air
                if (globalPos[0] > curHydraulicHead_) {
                    // try to model the initial water distribution
                    // above the hydraulic head due to the capillary
                    // pressure
#ifdef USE_NODE_PARAMETERS
                    dest[snIdx] = 1.0 - ParkerLenhard::Sw(ParentType::vertexState(0), - pH);
#else
                    const ElementState &cs = ParentType::elementState(element);
                    dest[snIdx] = 1.0 - ParkerLenhard::Sw(cs, - pH);
#endif

                    dest[snIdx] = std::min((Scalar) 1.0, dest[snIdx]);
                    dest[snIdx] = std::max((Scalar) 0.0, dest[snIdx]);

#ifdef USE_NODE_PARAMETERS
                    dest[pWIdx] = -ParkerLenhard::pC(ParentType::vertexState(0),
                                                     1 - dest[snIdx]);
#else
                    dest[pWIdx] = -ParkerLenhard::pC(cs,
                                                     1 - dest[snIdx]);
#endif


                }
                else {
                    dest[pWIdx] = pH;
                    dest[snIdx] = 0.0;
                }
            }

    public:
        // Returns the type of an boundary condition for the wetting
        // phase pressure at a element face
        void boundaryTypes(BoundaryTypeVector         &values,
                           const Element              &element,
                           const FVElementGeometry    &fvElemGeom,
                           const IntersectionIterator &isIt,
                           int                         scvIdx,
                           int                         boundaryFaceIdx) const
            {
//                const GlobalPosition &globalPos
//                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
//                const LocalPosition &localPos
//                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

                values[pWIdx] = Dune::BoundaryConditions::dirichlet;
                values[snIdx] = Dune::BoundaryConditions::dirichlet;
            }

        //! Evaluate a neumann boundary condition
        void neumann(SolutionVector             &values,
                     const Element              &element,
                     const FVElementGeometry    &fvElemGeom,
                     const IntersectionIterator &isIt,
                     int                         scvIdx,
                     int                         boundaryFaceIdx) const
            {
//                const GlobalPosition &globalPos
//                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
//                const LocalPosition &localPos
//                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

                values[pWIdx] = 0;
                values[snIdx] = 0;
            }


        //! Evaluate a dirichlet boundary condition at a element vert
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


#if defined USE_NODE_PARAMETERS
                if (onUpperBoundary(globalPos)) {
                    values[pWIdx] = -maxPc_;
                    values[snIdx] = 1;
                }
                else { // onLowerBoundary(pos)
                    Scalar h = curHydraulicHead_ - globalPos[0];
                    Scalar pH = hydrostaticPressure_(h);
                    values[pWIdx] = pH;
                    values[snIdx] = 0;
                }
#else
                initial_(dest, globalPos);
#endif
            }

        //! evaluate the mass injection rate of the fluids for a BOX
        //! sub control volume
        void source(SolutionVector          &dest,
                    const Element           &element,
                    const FVElementGeometry &dualElement,
                    int                      scvIdx)
            {
                dest[pWIdx] = dest[snIdx] = 0;
            }
        ///////////////////////////////////
        // End of problem specific stuff
        ///////////////////////////////////

        ///////////////////////////////////
        // Strings pulled by the Lenhard newton controller
        ///////////////////////////////////
        //! called by the LenhardNewtonContoller when the newton method
        //! is started.
        void newtonBegin()
            {
#if defined LENHARD_WRITE_NEWTON_STEPS
                convergenceWriter_ =
                    new VtkMultiWriter((boost::format("lenhard-convergence-t=%.2f-dt=%.2f")
                                        %timeManager_.time()%timeManager_.stepSize()).str());
#endif // LENHARD_WRITE_NEWTON_STEPS
            }

        //! called by the LenhardNewtonContoller when a newton step is
        //! finished.
        void newtonEndStep(SpatialFunction &u, SpatialFunction &uOld)
            {
#if defined LENHARD_WRITE_NEWTON_STEPS
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
#endif // LENHARD_WRITE_NEWTON_STEPS
            }

        //! called by the LenhardNewtonContoller when the newton method
        //! is finished.
        void newtonEnd()
            {
#if defined LENHARD_WRITE_NEWTON_STEPS
                delete convergenceWriter_;
#endif // LENHARD_WRITE_NEWTON_STEPS
            }
        ///////////////////////////////////
        // end of LenhardNewtonController stuff
        ///////////////////////////////////

    private:
        //! called when the solution for a time step has been
        //! computed.
#if LENHARD_EXPERIMENT == 1
        void updateEpisode_()
            {
                int epiIdx = timeManager_.episodeIndex();

                int i = 0;
                int k = 12;
                if (epiIdx == i) {
                    // initial episode, we start at t=-3hours
                    curHydraulicHead_ = 0.67;
                    timeManager_.startNextEpisode(-3*60*60);
                    timeManager_.setStepSize(initialTimeStepSize_);
                    return;
                }
                ++ i;

                if (timeManager_.time() >= endTime_) {
                    timeManager_.setFinished();
                    return;
                }
                else if (!timeManager_.episodeIsOver())
                    return;

                if (epiIdx < i + k) {
                    // reduce the hydraulic head by 5cm and wait
                    // for 10 minutes
                    curHydraulicHead_ -= 0.05;
                    timeManager_.startNextEpisode(10*60);
                    timeManager_.setStepSize(1);
                    return;
                };
                i += k;

                if (epiIdx == i) {
                    // wait until we reach simulation time t=53 hours
                    timeManager_.startNextEpisode(53*60*60 - timeManager_.time());
                    timeManager_.setStepSize(1);
                    return;
                }
                ++i;

                if (epiIdx < i + k) {
                    // increase the hydraulic head by 5cm and wait
                    // for 10 minutes
                    curHydraulicHead_ += 0.05;
                    timeManager_.startNextEpisode(10*60);
                    timeManager_.setStepSize(1);
                    return;
                }
                i += k;

                if (epiIdx == i) {
                    timeManager_.startNextEpisode(20*60*60);
                    timeManager_.setStepSize(1);
                    return;
                }

                timeManager_.setFinished();
            }

#elif LENHARD_EXPERIMENT == 2
        void updateEpisode_()
            {
                int epiIdx = timeManager_.episodeIndex();

                int i = 0;
                if (epiIdx == i) {
                    // initial episode, we start at t=-3hours
                    curHydraulicHead_ = 0.72;
                    timeManager_.setTime(-10*60*60, 0);
                    timeManager_.startNextEpisode(10*60*60);
                    timeManager_.setStepSize(initialTimeStepSize_);
                    return;
                }
                ++ i;

                if (timeManager_.time() >= endTime_) {
                    timeManager_.setFinished();
                    return;
                }
                else if (!timeManager_.episodeIsOver())
                    return;

                int k = 13;
                if (epiIdx < i + k) {
                    // reduce the hydraulic head by 5cm and wait
                    // for 10 minutes
                    curHydraulicHead_ -= 0.05;
                    timeManager_.startNextEpisode(10*60);
                    timeManager_.setStepSize(10.0);
                    return;
                };
                i += k;

                if (epiIdx == i) {
                    // wait until we reach simulation time t=3 hours
                    timeManager_.startNextEpisode(3*60*60 - timeManager_.time());
                    timeManager_.setStepSize(10.0);
                    return;
                }
                ++i;

                k = 7;
                if (epiIdx < i + k) {
                    // increase the hydraulic head by 5cm and wait
                    // for 10 minutes
                    curHydraulicHead_ += 0.05;
                    timeManager_.startNextEpisode(10*60);
                    timeManager_.setStepSize(10.0);
                    return;
                }
                i += k;

                if (epiIdx == i ) {
                    // wait until we reach simulation time t=3 hours
                    timeManager_.startNextEpisode(5*60*60 - timeManager_.time());
                    timeManager_.setStepSize(10.0);
                    return;
                }
                ++i;

                k = 5;
                if (epiIdx < i + k) {
                    // decrease the hydraulic head by 5cm and wait
                    // for 10 minutes
                    curHydraulicHead_ -= 0.05;
                    timeManager_.startNextEpisode(10*60);
                    timeManager_.setStepSize(10.0);
                    return;
                }
                i += k;

                if (epiIdx == i ) {
                    // wait until we reach simulation time t=6.67 hours
                    timeManager_.startNextEpisode(6.67*60*60 - timeManager_.time());
                    timeManager_.setStepSize(10.0);
                    return;
                }
                ++i;

                k = 11;
                if (epiIdx < i + k) {
                    // increase the hydraulic head by 5cm and wait
                    // for 10 minutes
                    curHydraulicHead_ += 0.05;
                    timeManager_.startNextEpisode(10*60);
                    timeManager_.setStepSize(60);
                    return;
                }
                i += k;

                if (epiIdx == i ) {
                    // wait until we reach simulation time t=10 hours
                    timeManager_.startNextEpisode(10*60*60 - timeManager_.time());
                    timeManager_.setStepSize(10.0);
                    return;
                }
                ++i;

                timeManager_.setFinished();
            }
#endif

        // write results to the output files
        void writeCurrentResult_()
            {
                // write some points to stdout so that they can be
                // later filtered and used for plotting them into a
                // SVG
                writeCsv_();

                // write the actual result into a VTK dataset
                resultWriter_.beginTimestep(timeManager_.time(),
                                            ParentType::grid().leafView());

                resultWriter_.addScalarVertexFunction("Sn",
                                                      model_.currentSolution(),
                                                      snIdx);
                resultWriter_.addScalarVertexFunction("Pw",
                                                      model_.currentSolution(),
                                                      pWIdx);

                SpatialFunction globResidual(ParentType::grid());
                model_.evalGlobalResidual(globResidual);
                resultWriter_.addScalarVertexFunction("global defect Sn",
                                                      globResidual,
                                                      snIdx);
                resultWriter_.addScalarVertexFunction("global defect Pw",
                                                      globResidual,
                                                      pWIdx);

                resultWriter_.endTimestep();
            }

#if defined LENHARD_WRITE_NEWTON_STEPS
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

                convergenceWriter_->addScalarVertexFunction("Sn",
                                                            u,
                                                            ParentType::vertMap(),
                                                            snIdx);
                convergenceWriter_->addScalarVertexFunction("Pw",
                                                            u,
                                                            ParentType::vertMap(),
                                                            pWIdx);
                convergenceWriter_->addScalarVertexFunction("difference Sn",
                                                            diff,
                                                            ParentType::vertMap(),
                                                            snIdx);
                convergenceWriter_->addScalarVertexFunction("difference Pw",
                                                            diff,
                                                            ParentType::vertMap(),
                                                            pWIdx);

/*                writeVertexFields_(*convergenceWriter_, u);
                writeElementFields_(*convergenceWriter_, u);
*/
            };
#endif // LENHARD_WRITE_NEWTON_STEPS

        // called whenever a solution for a timestep has been computed.
        void updateDomain_()
            {
#ifdef USE_HYSTERESIS
#ifdef USE_NODE_PARAMETERS
                VertexIterator it = ParentType::vertexBegin();
                VertexIterator endit = ParentType::vertexEnd();
                for (; it != endit; ++it) {
                    int vertIdx = ParentType::vertexIdx(*it);
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
                    Scalar Sn = model_.currentSolution().evallocal(snIdx,
                                                     *it,
                                                     localPos);
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(elementState(elementIdx(*it)), Sw);
                }
#endif
#endif
            }

        void writeCsv_()
            {
                // write out the points at 67,57,47,37 and 27 cm in
                // order to allow them being plotted
                ElementIterator it = ParentType::elementBegin();
                ElementIterator endIt = ParentType::elementEnd();
#if LENHARD_EXPERIMENT == 1
                Scalar points[] = { 0.27, 0.37, 0.47, 0.57, 0.67, 1e100 };
#elif LENHARD_EXPERIMENT == 2
                Scalar points[] = { 0.30, 0.40, 0.50, 0.60, 0.70, 1e100 };
#endif
                int curI = 0;
                for (; it != endIt; ++it) {
                    LocalPosition localPos = it->geometry().local(points[curI]);
                    if (0 <= localPos[0] && localPos[0] < 1.0) {
                        // evaluate the solution for the non-wetting
                        // saturation at this point
                        Scalar Sn = model_.currentSolution().evallocal(snIdx,
                                                         *it,
                                                         localPos);

                        // print the water saturation vs time at the
                        // given point in order to allow it to be
                        // plotted
                        std::cout << boost::format("snip: pos=%.02f\n")%points[curI];
                        std::cout << timeManager_.time()/3600 << " " << 1.0 - Sn << "\n";
                        std::cout << "snap\n";

                        // print the capillary pressure vs water
                        // saturtion at some vert of the element which
                        // contains the current point of interest
                        const VertexState &vs = ParentType::vertexState(ParentType::elementIdx(*it));
                        localPos[0] = 0;
                        Sn = model_.currentSolution().evallocal(snIdx,
                                                                *it,
                                                                localPos);
                        std::cout << boost::format("snip: pC_pos=%.02f\n")%points[curI];
                        std::cout << 1.0 - Sn << " " << ParkerLenhard::pC(vs, 1.0 - Sn) << "\n";
                        std::cout << "snap\n";

                        ++curI;
                    }
                }
            }

        Scalar hydrostaticPressure_(Scalar hydraulicHead) const
            {
                return -ParentType::densityW()*
                        ParentType::gravity()[0]*
                        hydraulicHead;
            };


        // simulated time control stuff
        TimeManager     timeManager_;
        TimeIntegration timeIntegration_;
        Scalar          initialTimeStepSize_;
        Scalar          endTime_;

        Scalar curHydraulicHead_;
        Scalar maxPc_;

        Model            model_;
        NewtonMethod     newtonMethod_;
        NewtonController newtonCtl_;
        VtkMultiWriter   resultWriter_;

#if defined LENHARD_WRITE_NEWTON_STEPS
        VtkMultiWriter *convergenceWriter_;
#endif // LENHARD_WRITE_NEWTON_STEPS
    };
}
}


#endif
