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

#include <dumux/new_models/pwsn/pwsnboxmodel.hh>
#include <dumux/new_models/pwsn/pwsnnewtoncontroller.hh>

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
            PrimaryVariables = BoxTraits::PrimaryVariables,
            PwIndex = PwSnTraits::PwIndex,
            SnIndex = PwSnTraits::SnIndex
        };
        
        // copy some types from the traits for convenience
        typedef typename DomainTraits::Scalar                     Scalar;
        typedef typename DomainTraits::Grid                       Grid;
        typedef typename DomainTraits::Cell                       Cell;
        typedef typename DomainTraits::CellIterator               CellIterator;
        typedef typename DomainTraits::CellReferenceElement       CellReferenceElement;
        typedef typename DomainTraits::CellReferenceElements      CellReferenceElements;
        typedef typename DomainTraits::Node                       Node;
        typedef typename DomainTraits::NodeIterator               NodeIterator;
        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter IntersectionIteratorGetter;
        typedef typename DomainTraits::LocalCoord                 LocalCoord;
        typedef typename DomainTraits::WorldCoord                 WorldCoord;

        typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction               SpatialFunction;
        typedef typename BoxTraits::UnknownsVector                UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

        typedef typename MaterialTraits::ParkerLenhard            ParkerLenhard;
        typedef typename MaterialTraits::CellState                CellState;
        typedef typename MaterialTraits::NodeState                NodeState;

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
            : _model(*this),
              _newtonMethod(_model),
              _newtonCtl(*this),
#if LENHARD_EXPERIMENT == 1
              _resultWriter("LenhardExp1")
#elif LENHARD_EXPERIMENT == 2
              _resultWriter("LenhardExp2")
#endif
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

#ifdef USE_NODE_PARAMETERS
                _maxPc = ParkerLenhard::pC(ParentType::nodeState(0), 
                                           0.0);
#else // !USE_NODE_PARAMETERS
                _maxPc = ParkerLenhard::pC(ParentType::cellState(0), 
                                           0.0);
#endif
                _initialTimeStepSize = initialTimeStepSize;
                _endTime = endTime;
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

                _timeManager.runSimulation(*this);

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
                _updateEpisode();

                 // set the initial solution
                _model.initial();

                // write the inital solution to disk
                _writeCurrentResult();
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
                // integration modifies the _curTimeStepSize if it
                // thinks that's appropriate. (IMHO, this is an
                // incorrect abstraction, the simulation
                // controller is responsible for adapting the time
                // step size!)
                _timeIntegration.execute(*this,
                                         _timeManager.time(),
                                         stepSize,
                                         nextStepSize,
                                         1e100, // firstDt or maxDt, TODO: WTF?
                                         _endTime,
                                         1.0); // CFL factor (not relevant since we use implicit euler)
            };


        //! called by the TimeIntegration::execute function to let
        //! time pass.
        void updateModel(Scalar &dt, Scalar &nextDt)
            {
                _model.update(dt, nextDt, _newtonMethod, _newtonCtl);
            }

        //! called by the TimeManager whenever a solution for a
        //! timestep has been computed
        void timestepDone()
            {
                // write the current result to disk
                _writeCurrentResult();

                // update the domain with the current solution
                _updateDomain();

                // change the episode of the simulation if necessary
                _updateEpisode();
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
            { return _timeManager.stepSize(); }

        //! Set the time step size in seconds.
        void setTimeStepSize(Scalar dt) 
            { return _timeManager.setStepSize(dt); }

        //! evaluate the initial condition for a node
        void initial(UnknownsVector &dest,
                     const Cell &cell,
                     WorldCoord pos,
                     LocalCoord posLocal)
            {
                Scalar h = _curHydraulicHead - pos[0];
                Scalar pH = _hydrostaticPressure(h);

                // initially the lower 67cm are filled with water the
                // remaining 5 with air
                if (pos[0] > _curHydraulicHead) {
                    // try to model the initial water distribution
                    // above the hydraulic head due to the capillary
                    // pressure
#ifdef USE_NODE_PARAMETERS
                    dest[SnIndex] = 1.0 - ParkerLenhard::Sw(ParentType::nodeState(0), - pH);
#else
                    const CellState &cs = ParentType::cellState(cell);
                    dest[SnIndex] = 1.0 - ParkerLenhard::Sw(cs, - pH);
#endif
                    
                    dest[SnIndex] = std::min((Scalar) 1.0, dest[SnIndex]);
                    dest[SnIndex] = std::max((Scalar) 0.0, dest[SnIndex]);
                    
#ifdef USE_NODE_PARAMETERS
                    dest[PwIndex] = -ParkerLenhard::pC(ParentType::nodeState(0),
                                                       1 - dest[SnIndex]);
#else
                    dest[PwIndex] = -ParkerLenhard::pC(cs,
                                                       1 - dest[SnIndex]);
#endif
                    

                }
                else {
                    dest[PwIndex] = pH;
                    dest[SnIndex] = 0.0;
                }
            }


        // Returns the type of an boundary contition for the wetting
        // phase pressure at a cell face
        void boundaryTypes(BoundaryTypeVector &dest,
                           const Cell &cell,
                           const IntersectionIterator &face,
                           const WorldCoord &pos,
                           const LocalCoord &localPos)

            {
                dest[PwIndex] = Dune::BoundaryConditions::dirichlet;
                dest[SnIndex] = Dune::BoundaryConditions::dirichlet;
            }

        //! Evaluate a neumann boundary condition
        void neumann(UnknownsVector &dest,
                     const Cell &cell,
                     const IntersectionIterator &face,
                     const WorldCoord &pos,
                     const LocalCoord &localPos)
            {
                dest[PwIndex] = 0;
                dest[SnIndex] = 0;
            }


        //! Evaluate a dirichlet boundary condition at a cell node
        void dirichlet(UnknownsVector &dest,
                       const Cell &cell,
                       int nodeIdx,
                       int globalIdx)
            {
                const LocalCoord &localPos = cell.geometry()[nodeIdx];
                WorldCoord pos = cell.geometry().global(localPos);
                
#if defined USE_NODE_PARAMETERS
                if (onUpperBoundary(pos)) {
                    dest[PwIndex] = -_maxPc;
                    dest[SnIndex] = 1;
                }
                else { // onLowerBoundary(pos) 
                    Scalar h = _curHydraulicHead - pos[0];
                    Scalar pH = _hydrostaticPressure(h);
                    dest[PwIndex] = pH;
                    dest[SnIndex] = 0;
                }
#else
                initial(dest, cell, pos, localPos);
#endif
            }

        //! evaluate the mass injection rate of the fluids for a BOX
        //! sub control volume
        void sourceTerm(UnknownsVector &dest,
                        const Cell &cell,
                        const FVElementGeometry &dualCell,
                        int subControlVolumeId)
            {
                dest[PwIndex] = dest[SnIndex] = 0;
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
                _convergenceWriter =
                    new VtkMultiWriter((boost::format("lenhard-convergence-t=%.2f-dt=%.2f")
                                        %_timeManager.time()%_timeManager.stepSize()).str());
#endif // LENHARD_WRITE_NEWTON_STEPS
            }

        //! called by the LenhardNewtonContoller when a newton step is
        //! finished.
        void newtonEndStep(SpatialFunction &u, SpatialFunction &uOld)
            {
#if defined LENHARD_WRITE_NEWTON_STEPS
                if (_newtonCtl.newtonNumSteps() == 1) {
                    _convergenceWriter->beginTimestep(0,
                                                      ParentType::grid().leafView());
                    _writeConvergenceFields(uOld, uOld);
                    _convergenceWriter->endTimestep();
                }


                _convergenceWriter->beginTimestep(_newtonCtl.newtonNumSteps(),
                                                  ParentType::grid().leafView());
                _writeConvergenceFields(u, uOld);
                _convergenceWriter->endTimestep();
#endif // LENHARD_WRITE_NEWTON_STEPS
            }

        //! called by the LenhardNewtonContoller when the newton method
        //! is finished.
        void newtonEnd()
            {
#if defined LENHARD_WRITE_NEWTON_STEPS
                delete _convergenceWriter;
#endif // LENHARD_WRITE_NEWTON_STEPS
            }
        ///////////////////////////////////
        // end of LenhardNewtonController stuff
        ///////////////////////////////////

    private:
        //! called when the solution for a time step has been
        //! computed.
#if LENHARD_EXPERIMENT == 1
        void _updateEpisode()
            {
                int epiIndex = _timeManager.episodeIndex();

                int i = 0;
                int k = 12;
                if (epiIndex == i) {
                    // initial episode, we start at t=-3hours
                    _curHydraulicHead = 0.67;
                    _timeManager.startNextEpisode(-3*60*60);
                    _timeManager.setStepSize(_initialTimeStepSize);
                    return;
                }
                ++ i;

                if (_timeManager.time() >= _endTime) {
                    _timeManager.setFinished();
                    return;
                }
                else if (!_timeManager.episodeIsOver())
                    return;

                if (epiIndex < i + k) {
                    // reduce the hydraulic head by 5cm and wait
                    // for 10 minutes
                    _curHydraulicHead -= 0.05;
                    _timeManager.startNextEpisode(10*60);
                    _timeManager.setStepSize(1);
                    return;
                };
                i += k;
                    
                if (epiIndex == i) {
                    // wait until we reach simulation time t=53 hours
                    _timeManager.startNextEpisode(53*60*60 - _timeManager.time());
                    _timeManager.setStepSize(1);
                    return;
                }
                ++i;
                
                if (epiIndex < i + k) {
                    // increase the hydraulic head by 5cm and wait
                    // for 10 minutes
                    _curHydraulicHead += 0.05;
                    _timeManager.startNextEpisode(10*60);
                    _timeManager.setStepSize(1);
                    return;
                }
                i += k;
                
                if (epiIndex == i) {
                    _timeManager.startNextEpisode(20*60*60);
                    _timeManager.setStepSize(1);
                    return;
                }
                
                _timeManager.setFinished();
            }

#elif LENHARD_EXPERIMENT == 2
        void _updateEpisode()
            {
                int epiIndex = _timeManager.episodeIndex();

                int i = 0;
                if (epiIndex == i) {
                    // initial episode, we start at t=-3hours
                    _curHydraulicHead = 0.72;
                    _timeManager.setTime(-10*60*60, 0);
                    _timeManager.startNextEpisode(10*60*60);
                    _timeManager.setStepSize(_initialTimeStepSize);
                    return;
                }
                ++ i;

                if (_timeManager.time() >= _endTime) {
                    _timeManager.setFinished();
                    return;
                }
                else if (!_timeManager.episodeIsOver())
                    return;

                int k = 13;
                if (epiIndex < i + k) {
                    // reduce the hydraulic head by 5cm and wait
                    // for 10 minutes
                    _curHydraulicHead -= 0.05;
                    _timeManager.startNextEpisode(10*60);
                    _timeManager.setStepSize(10.0);
                    return;
                };
                i += k;
                    
                if (epiIndex == i) {
                    // wait until we reach simulation time t=3 hours
                    _timeManager.startNextEpisode(3*60*60 - _timeManager.time());
                    _timeManager.setStepSize(10.0);
                    return;
                }
                ++i;
                
                k = 7;
                if (epiIndex < i + k) {
                    // increase the hydraulic head by 5cm and wait
                    // for 10 minutes
                    _curHydraulicHead += 0.05;
                    _timeManager.startNextEpisode(10*60);
                    _timeManager.setStepSize(10.0);
                    return;
                }
                i += k;
                
                if (epiIndex == i ) {
                    // wait until we reach simulation time t=3 hours
                    _timeManager.startNextEpisode(5*60*60 - _timeManager.time());
                    _timeManager.setStepSize(10.0);
                    return;
                }
                ++i;

                k = 5;
                if (epiIndex < i + k) {
                    // decrease the hydraulic head by 5cm and wait
                    // for 10 minutes
                    _curHydraulicHead -= 0.05;
                    _timeManager.startNextEpisode(10*60);
                    _timeManager.setStepSize(10.0);
                    return;
                }
                i += k;

                if (epiIndex == i ) {
                    // wait until we reach simulation time t=6.67 hours
                    _timeManager.startNextEpisode(6.67*60*60 - _timeManager.time());
                    _timeManager.setStepSize(10.0);
                    return;
                }
                ++i;

                k = 11;
                if (epiIndex < i + k) {
                    // increase the hydraulic head by 5cm and wait
                    // for 10 minutes
                    _curHydraulicHead += 0.05;
                    _timeManager.startNextEpisode(10*60);
                    _timeManager.setStepSize(60);
                    return;
                }
                i += k;
                
                if (epiIndex == i ) {
                    // wait until we reach simulation time t=10 hours
                    _timeManager.startNextEpisode(10*60*60 - _timeManager.time());
                    _timeManager.setStepSize(10.0);
                    return;
                }
                ++i;

                _timeManager.setFinished();
            }
#endif

        // write results to the output files
        void _writeCurrentResult()
            {
                // write some points to stdout so that they can be
                // later filtered and used for plotting them into a
                // SVG
                _writeCsv();
                
                // write the actual result into a VTK dataset
                _resultWriter.beginTimestep(_timeManager.time(),
                                            ParentType::grid().leafView());

                _resultWriter.addScalarVertexFunction("Sn",
                                                      _model.currentSolution(),
                                                      SnIndex);
                _resultWriter.addScalarVertexFunction("Pw",
                                                      _model.currentSolution(),
                                                      PwIndex);

                SpatialFunction globResidual(ParentType::grid());
                _model.evalGlobalResidual(globResidual);
                _resultWriter.addScalarVertexFunction("global defect Sn",
                                                      globResidual,
                                                      SnIndex);
                _resultWriter.addScalarVertexFunction("global defect Pw",
                                                      globResidual,
                                                      PwIndex);

                _resultWriter.endTimestep();
            }

#if defined LENHARD_WRITE_NEWTON_STEPS
        void _writeConvergenceFields(SpatialFunction &u, SpatialFunction &uOld)
            {
                SpatialFunction diff(ParentType::grid());
                for (int i=0; i < (*diff).size(); ++i) {
                    (*diff)[i] = (*u)[i] - (*uOld)[i];
                    if ((*diff)[i][0] > 1e12)
                        (*diff)[i][0] = 1e12;
                    if ((*diff)[i][1] > 1e12)
                        (*diff)[i][1] = 1e12;
                }
                
                _convergenceWriter->addScalarVertexFunction("Sn",
                                                            u,
                                                            ParentType::nodeMap(),
                                                            SnIndex);
                _convergenceWriter->addScalarVertexFunction("Pw",
                                                            u,
                                                            ParentType::nodeMap(),
                                                            PwIndex);
                _convergenceWriter->addScalarVertexFunction("difference Sn",
                                                            diff,
                                                            ParentType::nodeMap(),
                                                            SnIndex);
                _convergenceWriter->addScalarVertexFunction("difference Pw",
                                                            diff,
                                                            ParentType::nodeMap(),
                                                            PwIndex);
                
/*                _writeVertexFields(*_convergenceWriter, u);
                _writeCellFields(*_convergenceWriter, u);
*/
            };
#endif // LENHARD_WRITE_NEWTON_STEPS

        // called whenever a solution for a timestep has been computed.
        void _updateDomain()
            {
#ifdef USE_HYSTERESIS
#ifdef USE_NODE_PARAMETERS
                NodeIterator it = ParentType::nodeBegin();
                NodeIterator endit = ParentType::nodeEnd();
                for (; it != endit; ++it) {
                    int vertIdx = ParentType::nodeIndex(*it);
                    Scalar Sn = (*_model.currentSolution())[vertIdx][SnIndex];
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(nodeState(*it), Sw);
                }
#else
                // update the parker-lenhard hystersis model
                // for all cells
                CellIterator it = ParentType::grid().template leafbegin<0>();
                CellIterator endit = ParentType::grid().template leafend<0>();
                for (; it != endit; ++it) {
                    // get the barycenter of the current cell
                    const Dune::GeometryType &geoType = it->type();
                    const LocalCoord &localPos =
                        CellReferenceElements::general(geoType).position(0,0);

                    // evaluate the solution for the non-wetting
                    // saturation at this point
                    Scalar Sn = _model.currentSolution().evallocal(SnIndex,
                                                     *it,
                                                     localPos);
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(cellState(cellIndex(*it)), Sw);
                }
#endif
#endif
            }
        
        void _writeCsv()
            {
                // write out the points at 67,57,47,37 and 27 cm in
                // order to allow them being plotted
                CellIterator it = ParentType::cellBegin();
                CellIterator endIt = ParentType::cellEnd();
#if LENHARD_EXPERIMENT == 1
                Scalar points[] = { 0.27, 0.37, 0.47, 0.57, 0.67, 1e100 };
#elif LENHARD_EXPERIMENT == 2
                Scalar points[] = { 0.30, 0.40, 0.50, 0.60, 0.70, 1e100 };
#endif
                int curI = 0;
                for (; it != endIt; ++it) {
                    LocalCoord localPos = it->geometry().local(points[curI]);
                    if (0 <= localPos[0] && localPos[0] < 1.0) {
                        // evaluate the solution for the non-wetting
                        // saturation at this point
                        Scalar Sn = _model.currentSolution().evallocal(SnIndex,
                                                         *it,
                                                         localPos);

                        // print the water saturation vs time at the
                        // given point in order to allow it to be
                        // plotted
                        std::cout << boost::format("snip: pos=%.02f\n")%points[curI];
                        std::cout << _timeManager.time()/3600 << " " << 1.0 - Sn << "\n"; 
                        std::cout << "snap\n";

                        // print the capillary pressure vs water
                        // saturtion at some node of the cell which
                        // contains the current point of interest
                        const NodeState &vs = ParentType::nodeState(ParentType::cellIndex(*it));
                        localPos[0] = 0;
                        Sn = _model.currentSolution().evallocal(SnIndex,
                                                                *it,
                                                                localPos);
                        std::cout << boost::format("snip: pC_pos=%.02f\n")%points[curI];
                        std::cout << 1.0 - Sn << " " << ParkerLenhard::pC(vs, 1.0 - Sn) << "\n"; 
                        std::cout << "snap\n";
                        
                        ++curI;
                    }
                }
            }

        Scalar _hydrostaticPressure(Scalar hydraulicHead) 
            {
                return -ParentType::densityW()*
                        ParentType::gravity()[0]*
                        hydraulicHead;
            };


        // simulated time control stuff
        TimeManager     _timeManager;
        TimeIntegration _timeIntegration;
        Scalar          _initialTimeStepSize;
        Scalar          _endTime;

        Scalar _curHydraulicHead;
        Scalar _maxPc;

        Model            _model;
        NewtonMethod     _newtonMethod;
        NewtonController _newtonCtl;
        VtkMultiWriter   _resultWriter;

#if defined LENHARD_WRITE_NEWTON_STEPS
        VtkMultiWriter *_convergenceWriter;
#endif // LENHARD_WRITE_NEWTON_STEPS
    };
}
}


#endif
