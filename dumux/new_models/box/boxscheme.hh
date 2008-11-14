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
#ifndef DUMUX_BOX_SCHEME_HH
#define DUMUX_BOX_SCHEME_HH

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
     * \brief The base class for the BOX hybrid finite element/finite volume discretization scheme
     */
    template<class Implementation, 
             class BoxTraitsT, 
             class ProblemT,
             class LocalJacobianT>
    class BoxScheme
    {
        // copy the relevant problem specfific types from the problem
        // controller class
        typedef BoxScheme<Implementation, 
                          BoxTraitsT,
                          ProblemT, 
                          LocalJacobianT> ThisType;
        typedef ProblemT                  Problem;

    public:
        /*!
         * \brief The traits of the BOX scheme.
         *
         * This includes the shape functions to be used, etc.
         */
        typedef BoxTraitsT                     BoxTraits;
        /*!
         * \brief The traits of the spatial domain (grid type, etc)
         */
        typedef typename Problem::DomainTraits DomainTraits;
        
        /*!
         *  \brief This structure is Required to use models based on the BOX
         *         scheme in conjunction with the newton Method.
         */
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

        typedef typename BoxTraits::JacobianAssembler          JacobianAssembler;
        typedef typename BoxTraits::SpatialFunction            SpatialFunction;
        typedef typename SpatialFunction::RepresentationType   BoxFnRep;
        typedef typename BoxTraits::LocalFunction              LocalFunction;

        typedef typename BoxTraits::ShapeFunctionSetContainer  ShapeFunctionSetContainer;
        typedef typename ShapeFunctionSetContainer::value_type ShapeFunctionSet;

        typedef typename BoxTraits::BoundaryTypeVector  UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector  BoundaryTypeVector;

        typedef LocalJacobianT                          LocalJacobian;

        // some constants
        enum {
            PrimaryVariables = BoxTraits::PrimaryVariables,

            GridDim     = DomainTraits::GridDim,
            WorldDim    = DomainTraits::WorldDim
        };
        
    public:
        BoxScheme(Problem &prob, LocalJacobian &localJac)
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

        /*!
         * \brief Reference to the current solution.
         */
        const SpatialFunction &currentSolution() const
            { return _uCur; }

        /*!
         * \brief Reference to the current solution.
         */
        SpatialFunction &currentSolution()
            { return _uCur; }

        /*!
         * \brief Reference to the right hand side.
         */
        SpatialFunction &rightHandSide()
            { return _f; }

        /*!
         * \brief Reference to solution of the previous time step.
         */
        SpatialFunction &previousSolution()
            { return _uPrev; }

        /*!
         * \brief Reference to solution of the previous time step.
         */
        const SpatialFunction &previousSolution() const
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
        LocalJacobian &localJacobian()
            { return _localJacobian; }

        /*!
         * \brief Same as localJacobian(), included to ease porting.
         */
        LocalJacobian &getLocalJacobian() DUNE_DEPRECATED
            { return _localJacobian; }

        /*!
         * \brief Reference to the grid of the spatial domain.
         */
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
                // here (previously it was in the newton solver where it
                // belongs even less)
                int numRetries = 0;
                while (true)
                {
                    bool converged = solver.execute(*this->_asImp(),
                                                    controller);
                    nextDt = controller.suggestTimeStepSize(dt);
                    
                    if (converged)
                        break;
                    
                    if (numRetries >= 10)
                        DUNE_THROW(Dune::MathError,
                                   "Newton solver didn't converge after 10 timestep divisions. dt=" << dt);
                    ++numRetries;
                    _problem.setTimeStepSize(nextDt);
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
         * The global deflection of the mass balance from zero.
         */
        void evalGlobalResidual(SpatialFunction &globResidual)
            {
                (*globResidual) = Scalar(0.0);

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
                    const int numVertices = cell.template count<GridDim>();
                    LocalFunction localResidual(numVertices);
                    LocalFunction localU(numVertices);
                    LocalFunction localOldU(numVertices);

                    _localJacobian.setCurrentCell(cell);
                    _localJacobian.restrictToCell(localU, currentSolution());
                    _localJacobian.restrictToCell(localOldU, previousSolution());
                    _localJacobian.setParams(cell, localU, localOldU);
                    _localJacobian.evalLocalResidual(localResidual);

                    // loop over the cell's vertices, map them to the
                    // corresponding grid's node ids and add the
                    // cell's local residual at a node the global
                    // residual at this node.
                    const ShapeFunctionSet &shapeFns = BoxTraits::shapeFunctions(cell.geometry().type(), 1);
                    for(int localId=0; localId < shapeFns.size(); localId++)
                    {
                        int globalId = _problem.nodeIndex(cell,
                                                            shapeFns[localId].entity());
                        (*globResidual)[globalId] += localResidual[localId];
                    }
                }
            }


    protected:
        void _applyInitialSolution(SpatialFunction &u)
            {
                // iterate through leaf grid an evaluate c0 at cell center
                CellIterator it     = _problem.grid().template leafbegin<0>();
                CellIterator eendit = _problem.grid().template leafend<0>();
                for (; it != eendit; ++it)
                {
                    // loop over all shape functions of the current cell
                    const Cell& cell = *it;
                    const ShapeFunctionSet &shapeFnSet = BoxTraits::shapeFunctions(cell.geometry().type(), 1);
                    for (int i = 0; i < shapeFnSet.size(); i++) {
                        // get the local and global coordinates of the
                        // shape function's center (i.e. the node
                        // where it is 1 for Lagrange functions)
                        const LocalCoord &localPos = shapeFnSet[i].position();
                        WorldCoord globalPos = it->geometry().global(localPos);

                        // translate the local index of the center of
                        // the current shape function to the global
                        // node id
                        int globalId = _problem.nodeIndex(cell, shapeFnSet[i].entity());

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
                    const typename ShapeFunctionSetContainer::value_type &shapeFnSet = BoxTraits::shapeFunctions(geoType, 1);

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

                            // translate local node id to a global one
                            int globalId = _problem.nodeIndex(cell,
                                                                shapeFnSet[i].entity());

                            // actually evaluate the boundary
                            // condition for the current
                            // cell+face+node combo
                            _problem.dirichlet((*u)[globalId],
                                               cell,
                                               faceIt,
                                               globalPos,
                                               localPos);                            
                        }
                    }
                }
            };

        Implementation *_asImp() 
            { return static_cast<Implementation*>(this); } 
        const Implementation *_asImp() const
            { return static_cast<const Implementation*>(this); } 

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
