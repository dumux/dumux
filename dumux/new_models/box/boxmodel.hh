/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
#ifndef DUMUX_BOX_MODEL_HH
#define DUMUX_BOX_MODEL_HH

#include <dumux/new_models/box/boxjacobian.hh>
#include <dumux/auxiliary/basicdomain.hh>
#include <dumux/nonlinear/new_newtonmethod.hh>
#include <dumux/auxiliary/apis.hh>

#include <dune/istl/operators.hh>
#include <dune/disc/operators/p1operator.hh>

#include <boost/format.hpp>

namespace Dune
{

    /*!
     * \brief The base class for the BOX hybrid finite element/finite volume discretization model
     */
    template<class BoxTraitsT, class ProblemT, class LocalJacobianT>
    class BoxModel
    {
        // copy the relevant problem specfific types from the problem
        // controller class
        typedef BoxModel<BoxTraitsT, ProblemT, LocalJacobianT> ThisType;
        typedef ProblemT                                       Problem;

    public:
        typedef BoxTraitsT                     BoxTraits;
        typedef typename Problem::DomainTraits DomainTraits;
        
        // required to use the model in conjunction with the newton
        // method
        struct NewtonTraits {
            typedef LocalJacobianT                         LocalJacobian;
            typedef typename BoxTraits::SpatialFunction    Function;
            typedef typename BoxTraits::JacobianAssembler  JacobianAssembler;
            typedef typename DomainTraits::Scalar          Scalar;
        };
        
    private:
        // copy the types from the traits for convenience
        typedef typename DomainTraits::Scalar                      Scalar;
        typedef typename DomainTraits::Grid                        Grid;
        typedef typename DomainTraits::Cell                        Cell;
        typedef typename DomainTraits::CellReferenceElement        CellReferenceElement;
        typedef typename DomainTraits::CellReferenceElements       CellReferenceElements;
        typedef typename DomainTraits::CellIterator                CellIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter  IntersectionIteratorGetter;
        typedef typename DomainTraits::IntersectionIterator        IntersectionIterator;
        typedef typename DomainTraits::CoordScalar                 CoordScalar;
        typedef typename DomainTraits::WorldCoord                  WorldCoord;
        typedef typename DomainTraits::LocalCoord                  LocalCoord;

        typedef typename BoxTraits::JacobianAssembler      JacobianAssembler;
        typedef typename BoxTraits::SpatialFunction        SpatialFunction;
        typedef typename SpatialFunction::RepresentationType   BoxFnRep;
        typedef typename BoxTraits::LocalFunction          LocalFunction;
        typedef typename BoxTraits::ShapeFnSets            ShapeFnSets;
        typedef typename BoxTraits::ShapeFnSet             ShapeFnSet;

        typedef typename BoxTraits::BoundaryTypeVector  UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector  BoundaryTypeVector;

        typedef LocalJacobianT                          LocalJacobian;

        // some constants
        enum {
            NumUnknowns = BoxTraits::NumUnknowns,

            GridDim     = DomainTraits::GridDim,
            WorldDim    = DomainTraits::WorldDim
        };
        
    public:
        BoxModel(Problem &prob, LocalJacobian &localJac)
            : _problem(prob),
              _uCur(prob.grid()),
              _uPrev(prob.grid()),
              _f(prob.grid()),
              _jacAsm(prob.grid()),
              _localJacobian(localJac)
            {
                Api::require<Api::BasicDomainTraits, 
                             typename Problem::DomainTraits>();
//                Api::require<Api::PwSnBoxDomain>(prob);
            }

        /*!
         * \brief Apply the initial conditions to the model.
         */
        void initial()
            {
                _applyInitialSolution(_uCur);
                _applyDirichletBoundaries(_uCur);
                
                *_uPrev = *_uCur;
            }

        // current solution
        const SpatialFunction &currentSolution() const
            { return _uCur; }

        // current solution
        SpatialFunction &currentSolution()
            { return _uCur; }

        // right hand side (?)
        SpatialFunction &f()
            { return _f; }

        // last timestep's solution
        SpatialFunction &uOldTimeStep()
            { return _uPrev; }

        const SpatialFunction &uOldTimeStep() const
            { return _uPrev; }

        /*!
         * \brief Returns the operator assembler for the global jacobian of
         *        the problem.
         */
        JacobianAssembler &jacobianAssembler()
            { return _jacAsm; }

        /*!
         * \brief Returns the local jacobian which calculates the local
         *        stiffness matrix for an arbitrary cell.
         * 
         * The local stiffness matrices of the cell are used by
         * the jacobian assembler to produce a global linerization of the
         * problem.
         */
        LocalJacobian &getLocalJacobian()
            { return _localJacobian; }

        const Grid &grid()
            { return _problem.grid(); }

        /*!
         * \brief Try to progress simulation to the next timestep.
         */
        template<class NewtonMethod, class NewtonController>
        void update(Scalar &dt, Scalar &nextDt, NewtonMethod &solver, NewtonController &controller)
            {
                _localJacobian.setCurrentSolution(&_uCur);
                _localJacobian.setOldSolution(&_uPrev);
                
                _applyDirichletBoundaries(_uCur);
                
                // TODO/FIXME: timestep control doesn't really belong
                // here (before it was in the newton solver where it
                // belongs even less)
                int numRetries = 0;
                while (true)
                {
                    _localJacobian.setDt(dt);
                    bool converged = solver.execute(*this, controller);
                    nextDt = controller.suggestTimeStepSize(dt);
                    
                    if (converged)
                        break;
                    
                    if (numRetries >= 10)
                        DUNE_THROW(Dune::MathError,
                                   "Newton solver didn't converge after 10 timestep divisions. dt=" << dt);
                    ++numRetries;
                    dt = nextDt;
                    std::cout << boost::format("Newton didn't converge. Retrying with timestep of %f\n")%dt;
                }

                // make the current solution the previous one. we copy
                // the whole representation here, because the current
                // solution is usually a much better approximation of
                // the next time step than the previous one. This
                // usually causes the newton solver to converge much
                // faster.
                *_uPrev = *_uCur;
            }


        /*!
         * \brief Calculate the global residual.
         * 
         * The global difference of the result when
         * using an approximate solution from the right hand side.
         */
        void evalGlobalResidual(SpatialFunction &globResidual)
            {
                (*globResidual)=0;

                // iterate through leaf grid
                CellIterator it     = _problem.grid().template leafbegin<0>();
                CellIterator eendit = _problem.grid().template leafend<0>();
                for (; it != eendit; ++it)
                {
                    // tell the local jacobian which cell it should
                    // consider and evaluate the local residual for the
                    // cell. in order to do this we first have to
                    // evaluate the cell's local solutions for the
                    // current and the last timestep.
                    const Cell& cell = *it;
                    LocalFunction localResidual;
                    LocalFunction localU;
                    LocalFunction localOldU;
                    _localJacobian.setCurrentCell(cell);
                    _localJacobian.evalLocal(localU, currentSolution());
                    _localJacobian.evalLocal(localOldU, uOldTimeStep());
                    _localJacobian.evalLocalResidual(localResidual,
                                                     localU,
                                                     localOldU);

                    // loop over the cell's vertices, map them to the
                    // corresponding grid's vertex ids and add the
                    // cell's local residual at a vertex the global
                    // residual at this vertex.
                    const ShapeFnSet &shapeFnSet = ShapeFnSets::general(cell.geometry().type(),
                                                                          1);
                    for(int localId=0; localId < shapeFnSet.size(); localId++)
                    {
                        int globalId = _problem.vertexIndex(cell,
                                                            shapeFnSet[localId].entity());
                        (*globResidual)[globalId] += localResidual.atSubContVol[localId];
                    }
                }
            }


    private:
        void _applyInitialSolution(SpatialFunction &u)
            {
                // iterate through leaf grid an evaluate c0 at cell center
                CellIterator it     = _problem.grid().template leafbegin<0>();
                CellIterator eendit = _problem.grid().template leafend<0>();
                for (; it != eendit; ++it)
                {
                    // loop over all shape functions of the current cell
                    const Cell& cell = *it;
                    const ShapeFnSet &shapeFnSet = ShapeFnSets::general(cell.geometry().type(), 1);
                    for (int i = 0; i < shapeFnSet.size(); i++) {
                        // get the local and global coordinates of the
                        // shape function's center (i.e. the vertex
                        // where it is 1 for Lagrange functions)
                        const LocalCoord &localPos = shapeFnSet[i].position();
                        WorldCoord globalPos = it->geometry().global(localPos);

                        // translate the local index of the center of
                        // the current shape function to the global
                        // vertex id
                        int globalId = _problem.vertexIndex(cell, shapeFnSet[i].entity());

                        // use the problem controller to actually do
                        // the dirty work of nailing down the initial
                        // solution.
                        _problem.initial((*u)[globalId],
                                         cell,
                                         globalPos,
                                         localPos);
                    }
                }
            };


        void _applyDirichletBoundaries(SpatialFunction &u)
            {
                // set Dirichlet boundary conditions of the grid's
                // outer boundaries
                CellIterator cellIt     = _problem.grid().template leafbegin<0>();
                CellIterator cellEndIt  = _problem.grid().template leafend<0>();
                for (; cellIt != cellEndIt; ++cellIt)
                {
                    if (!cellIt->hasBoundaryIntersections())
                        continue;
                    
                    // get the current cell and its set of shape
                    // functions
                    const Cell& cell = *cellIt;
                    Dune::GeometryType geoType = cell.geometry().type();
                    const ShapeFnSet &shapeFnSet = ShapeFnSets::general(geoType, 1);

                    // locally evaluate the cell's boundary condition types
                    _localJacobian.assembleBoundaryCondition(cell);
                    
                    // loop over all faces of the cell
                    const IntersectionIterator &faceEndIt = IntersectionIteratorGetter::end(cell);
                    IntersectionIterator       faceIt = IntersectionIteratorGetter::begin(cell);
                    for (; faceIt != faceEndIt;  ++faceIt) {
                        
                        // loop over all shape functions of the cell
                        for (int i = 0; i < shapeFnSet.size(); i++)
                        {
                            if (_localJacobian.bc(i)[0] != Dune::BoundaryConditions::dirichlet)
                                // we ought to evaluate dirichlet
                                // boundary conditions, not
                                // something else!
                                continue;

                            // get the shape function's center in
                            // local and global coordinates
                            const LocalCoord &localPos = shapeFnSet[i].position();
                            WorldCoord globalPos = cell.geometry().global(localPos);

                            // translate local vertex id to a global one
                            int globalId = _problem.vertexIndex(cell,
                                                                shapeFnSet[i].entity());

                            // actually evaluate the boundary
                            // condition for the current
                            // cell+face+vertex combo
                            _problem.dirichlet((*u)[globalId],
                                               cell,
                                               faceIt,
                                               globalPos,
                                               localPos);                            
                        }
                    }
                }
            };

        // the problem we want to solve. defines the constitutive
        // relations, material laws, etc.
        Problem     &_problem;

        // the solution we are looking for

        // cur is the current solution, prev the solution of the
        // previous time step
        SpatialFunction _uCur;
        SpatialFunction _uPrev;

        // the right hand side (?)
        SpatialFunction  _f;
        // Operator assembler. Linearizes the problem at a specific
        // position using the local jacobian (?)
        JacobianAssembler _jacAsm;
        // calculates the jacobian matrix at a given position
        LocalJacobian    &_localJacobian;
    };
}

#endif
