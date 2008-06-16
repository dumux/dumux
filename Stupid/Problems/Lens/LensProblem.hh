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
#ifndef STUPID_LENSPROBLEM_HH
#define STUPID_LENSPROBLEM_HH

#include "LensDomain.hh"
#include "LensNewtonController.hh"

#include <Stupid/Models/Box/PwSn/PwSnBoxModel.hh>
#include <Stupid/Material/ParkerLenhard.hh>
#include <Stupid/Material/RegularizedVanGenuchten.hh>
#include <Stupid/TimeDisc/ImplicitEulerStep.hh>
#include <Stupid/Auxilary/VtkMultiWriter.hh>
#include <Stupid/Auxilary/TimeManager.hh>

#include <iostream>

#define LENS_WRITE_NEWTON_STEPS 0

namespace Stupid
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
        // the traits of the BOX method
        typedef typename Model::BoxTraits           BoxTraits;
        // the traits for the material relations
        typedef typename ParentType::MaterialTraits MaterialTraits;

    private:
        // some constants from the traits for convenience
        enum {
            NumUnknowns = BoxTraits::NumUnknowns,
            PwIndex = BoxTraits::PwIndex,
            SnIndex = BoxTraits::SnIndex
        };
        
        // copy some types from the traits for convenience
        typedef typename DomainTraits::Scalar                     Scalar;
        typedef typename DomainTraits::Grid                       Grid;
        typedef typename DomainTraits::Cell                       Cell;
        typedef typename DomainTraits::CellIterator               CellIterator;
        typedef typename DomainTraits::CellReferenceElement       CellReferenceElement;
        typedef typename DomainTraits::CellReferenceElements      CellReferenceElements;
        typedef typename DomainTraits::Vertex                     Vertex;
        typedef typename DomainTraits::VertexIterator             VertexIterator;
        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter IntersectionIteratorGetter;
        typedef typename DomainTraits::LocalCoord                 LocalCoord;
        typedef typename DomainTraits::WorldCoord                 WorldCoord;

        typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
        typedef typename BoxTraits::BoxFunction                   BoxFunction;
        typedef typename BoxTraits::UnknownsVector                UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

        typedef typename MaterialTraits::ParkerLenhard            ParkerLenhard;

        // episode control stuff
        enum Episode {
            ImbibEpisode,  // an episode where imbibition of water takes place
            DrainEpisode,  // an episode where drainage of water takes place
            WaitEpisode    // an episode with neither drainage nor imbibition
        };

        typedef TimeManager<Episode>                TimeManager;
        typedef Stupid::ImplicitEulerStep<ThisType> TimeIntegration;
        typedef VtkMultiWriter<Grid>                VtkMultiWriter;

        typedef LensNewtonController<Model, ThisType>    NewtonController;

    public:
        PwSnLensProblem(Scalar initialTimeStepSize, Scalar endTime)
            : _model(*this),
              _newtonCtl(*this),
              _resultWriter("lens")
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                _initialTimeStepSize = initialTimeStepSize;
                _endTime = endTime;
            };

        ~PwSnLensProblem()
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
                    boost::format("LensSimulation took %.3f seconds\n")
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
                // set the initial condition
                _model.initial();

                // start with a drainage for 30 ksec
                _timeManager.startNextEpisode(DrainEpisode, 30e3);
                _timeManager.setStepSize(_initialTimeStepSize);

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
                _model.update(dt, nextDt, _newtonCtl);
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

        //! evaluate the initial condition for a vertex
        void initial(UnknownsVector &dest,
                     const Cell &cell,
                     WorldCoord pos,
                     LocalCoord posLocal)
            {
/*                WorldCoord pos;
                ParentType::vertexPosition(pos,
                                            ParentType::vertex(cell, localVertIdx));
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
                    Sn = _lensMedium->Snr();
                else
                    Sn = _outerMedium->Snr();
#else
                Sn = 0;
#endif

                dest[0] = pw; // pw
                dest[1] = Sn; // Sn
            }


        // Returns the type of an boundary contition for the wetting
        // phase pressure at a cell face
        void boundaryTypes(BoundaryTypeVector &dest,
                           const Cell &cell,
                           const IntersectionIterator &face,
                           const WorldCoord &pos,
                           const LocalCoord &localPos)

            {
                // get the integration point of the boundary face in
                // world coodinates
//        WorldCoord &pos = dualCell.boundaryFace[dcBFIndex].ipGlobal;

                if (ParentType::onLeftBoundary(pos) ||
                    ParentType::onRightBoundary(pos))
                {
                    dest[0] = dest[1] = Dune::BoundaryConditions::dirichlet;

#ifndef USE_ORIG_PROB
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
        void neumann(UnknownsVector &dest,
                     const Cell &cell,
                     const IntersectionIterator &face,
                     const WorldCoord &pos,
                     const LocalCoord &localPos)
            {
                dest[PwIndex] = 0;
                dest[SnIndex] = 0;

                // get the integration point of the boundary face in
                // world coodinates
//        WorldCoord &pos = dualCell.boundaryFace[dcBFIndex].ipGlobal;

                if (ParentType::onUpperBoundary(pos)) {
#ifdef USE_ORIG_PROB
                    Scalar relPosX = (ParentType::upperRight()[0] - pos[0])/ParentType::width();
                    if (0.5 < relPosX && relPosX < 2.0/3.0)
                    {
                        dest[SnIndex] = -0.04;
                    }
#endif
                }
#ifndef USE_ORIG_PROB
                else if (ParentType::onLowerBoundary(pos)) {
                    dest[SnIndex] = 0.0;
                    if (_timeManager.episode() == DrainEpisode)
                        // drain water
                        dest[SnIndex] = -0.04;
                    else if (_timeManager.episode() == ImbibEpisode)
                        // imbibition of water
                        dest[SnIndex] = 0.04;
                }
#endif

            }


        //! Evaluate a dirichlet boundary condition at a vertex within
        //! an cell's face
        void dirichlet(UnknownsVector &dest,
                       const Cell &cell,
                       const IntersectionIterator &face,
                       const WorldCoord &pos,
                       const LocalCoord &localPos)
            {
                Scalar a, b;

                // get the integration point of the boundary face in
                // world coodinates
//        WorldCoord &pos = dualCell.boundaryFace[dcBFIndex].ipGlobal;

                if (ParentType::onLeftBoundary(pos))
                {
                    a = -(1 + 0.5/ParentType::height());
                    b = -a*ParentType::upperRight()[1];
                }
                else {
                    a = -1;
                    b = ParentType::upperRight()[1];
                }

                dest[PwIndex] = -ParentType::densityW()*
                                  ParentType::gravity()[1]*
                                  (a*pos[1] + b);
#if 0
                dest[SnIndex] = ParentType::outerMedium().Snr();
#else
                dest[SnIndex] = 0;
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

        //! called by the LensNewtonContoller when the newton method
        //! is started.
        void newtonBegin()
            {
#if LENS_WRITE_NEWTON_STEPS
                _convergenceWriter =
                    new VtkMultiWriter((boost::format("lens-convergence-t=%.2f-dt=%.2f")
                                         %_curTime%_curTimeStepSize).str());
#endif // LENS_WRITE_NEWTON_STEPS
            }

        //! called by the LensNewtonContoller when a newton step is
        //! finished.
        void newtonEndStep(BoxFunction &u, BoxFunction &uOld)
            {
#if LENS_WRITE_NEWTON_STEPS
                if (_newtonCtl.newtonNumSteps() == 1) {
                    _convergenceWriter->beginTimestep(0,
                                                      grid());
                    _writeConvergenceFields(uOld, uOld);
                    _convergenceWriter->endTimestep();
                }


                _convergenceWriter->beginTimestep(_newtonCtl.newtonNumSteps(),
                                                  grid());
                _writeConvergenceFields(u, uOld);
                _convergenceWriter->endTimestep();
#endif // LENS_WRITE_NEWTON_STEPS
            }

        //! called by the LensNewtonContoller when the newton method
        //! is finished.
        void newtonEnd()
            {
#if LENS_WRITE_NEWTON_STEPS
                delete _convergenceWriter;
#endif // LENS_WRITE_NEWTON_STEPS
            }

    private:
        // write results to the output files
        void _writeCurrentResult()
            {
                _resultWriter.beginTimestep(_timeManager.time(),
                                            ParentType::grid());
                _writeVertexFields(_resultWriter, _model.u());
                _writeCellFields(_resultWriter, _model.u());
                _resultWriter.endTimestep();
            }

#if LENS_WRITE_NEWTON_STEPS
        void _writeConvergenceFields(BoxFunction &u, BoxFunction &uOld)
            {
                PwSnFunction diff(grid());
                for (int i=0; i < (*diff).size(); ++i) {
                    (*diff)[i] = (*u)[i] - (*uOld)[i];
                    if ((*diff)[i][0] > 1e12)
                        (*diff)[i][0] = 1e12;
                    if ((*diff)[i][1] > 1e12)
                        (*diff)[i][1] = 1e12;
                }

                _writeScalarVertexField(*_convergenceWriter,
                                        "reduction Sn",
                                        diff,
                                        SnIndex);
                _writeScalarVertexField(*_convergenceWriter,
                                        "reduction Pw",
                                        diff,
                                        PwIndex);
                _writeVertexFields(*_convergenceWriter, u);
                _writeCellFields(*_convergenceWriter, u);
            };
#endif // LENS_WRITE_NEWTON_STEPS


        // appends a cell centered capillary pressure field to the
        // multi writer's current timestep
        void _writeCellFields(VtkMultiWriter &writer, BoxFunction &u)
            {
                // TODO
                /*
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
                int nCells =  ParentType::numCells();
                ScalarField *pC = writer.template createField<Scalar, 1>(nCells);
                ScalarField *dpC_dSw = writer.template createField<Scalar, 1>(nCells);

                CellIterator it = ParentType::grid().template leafbegin<0>();
                CellIterator endIt = ParentType::grid().template leafend<0>();
                for (; it != endIt; ++it) {
                    // extract the current solution's Sn component
                    const CellReferenceElement &refElem =
                        CellReferenceElements::general(it->geometry().type());
                    Scalar Sn = u.evallocal(SnIndex,
                                             *it,
                                             refElem.position(0,0));

                    // calculate the capillary pressure
                    CellState &eState = ParentType::cellState(*it);
                    int eIndex = ParentType::cellIndex(*it);
                    Scalar Sw = 1 - Sn;
                    (*pC)[eIndex] = ParentType::pC(eState, Sw);
                    (*dpC_dSw)[eIndex] = ParentType::dpC_dSw(eState, Sw);
                }
                writer.addCellData(pC, "capillary pressure");
                writer.addCellData(dpC_dSw, "dpC/dSw");
                */
            }


        // write the fields current solution into an VTK output
        // file.
        void _writeVertexFields(VtkMultiWriter &writer, BoxFunction &u)
            {
                writer.addScalarVertexFunction("Sn",
                                               u,
                                               ParentType::vertexMap(),
                                               SnIndex);
                writer.addScalarVertexFunction("Pw",
                                               u,
                                               ParentType::vertexMap(),
                                               PwIndex);

                BoxFunction globDefect(ParentType::grid());
                _model.evalGlobalDefect(globDefect);
                writer.addScalarVertexFunction("global defect Sn",
                                               globDefect,
                                               ParentType::vertexMap(),
                                               SnIndex);
                writer.addScalarVertexFunction("global defect Pw",
                                               globDefect,
                                               ParentType::vertexMap(),
                                               PwIndex);
            }

        // called whenever a solution for a timestep has been computed.
        void _updateDomain()
            {
#ifdef USE_HYSTERESIS
#ifdef USE_VERTEX_PARAMETERS
                VertexIterator it = ParentType::vertexBegin();
                VertexIterator endit = ParentType::vertexEnd();
                for (; it != endit; ++it) {
                    int vertIdx = ParentType::vertexIndex(*it);
                    Scalar Sn = (*_model.u())[vertIdx][SnIndex];
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(vertexState(*it), Sw);
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
                    Scalar Sn = _model.u().evallocal(SnIndex,
                                                      *it,
                                                      localPos);
                    Scalar Sw = 1 - Sn;
                    ParkerLenhard::updateState(cellState(cellIndex(*it)), Sw);
                }
#endif
#endif
            }

        void _updateEpisode()
            {
                Scalar    len = _timeManager.episodeLength();
                Episode   epi = _timeManager.episode();
                int        epiIndex = _timeManager.episodeIndex();

                if (_timeManager.time() >= _endTime) {
                    _timeManager.setFinished();
                    return;
                }
                
#ifndef USE_ORIG_PROB               
                if (!_timeManager.episodeIsOver())
                    return;


                switch (epiIndex) {
                    case 1:
                        _timeManager.startNextEpisode(WaitEpisode, 50e3);
                        return;
                    case 2:
                        _timeManager.startNextEpisode(ImbibEpisode, 10e3);
                        return;
                    case 3:
                        _timeManager.startNextEpisode(WaitEpisode, 50e3);
                        return;
                    case 4:
                        _timeManager.startNextEpisode(DrainEpisode, 20e3);
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

                _timeManager.startNextEpisode(epi, len);
#endif
            }

        // simulated time control stuff
        TimeManager     _timeManager;
        TimeIntegration _timeIntegration;
        Scalar          _initialTimeStepSize;
        Scalar          _endTime;

        Model            _model;
        NewtonController _newtonCtl;
        VtkMultiWriter   _resultWriter;
#if LENS_WRITE_NEWTON_STEPS
        VtkMultiWriter *_convergenceWriter;
#endif // LENS_WRITE_NEWTON_STEPS

    };
}
}


#endif
