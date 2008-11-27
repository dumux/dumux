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

#include <dumux/new_models/boxscheme/boxjacobian.hh>
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
        typedef typename DomainTraits::ReferenceElement            ReferenceElement;
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

        typedef typename BoxTraits::UnknownsVector      UnknownsVector;
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
            : problem_(prob),
              uCur_(prob.grid()),
              uPrev_(prob.grid()),
              f_(prob.grid()),
              jacAsm_(prob.grid()),
              localJacobian_(localJac)
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
                // initialize the static node data of the box jacobian
                this->localJacobian().setCurSolution(&uCur_);
                this->localJacobian().setOldSolution(&uPrev_);
                
                this->localJacobian().initStaticData();
                
                applyInitialSolution_(uCur_);
                applyDirichletBoundaries_(uCur_);

                *uPrev_ = *uCur_;

                // update the static node data with the initial solution
                this->localJacobian().updateStaticData(uCur_, uPrev_);              
            }

        /*!
         * \brief Reference to the current solution.
         */
        const SpatialFunction &currentSolution() const
            { return uCur_; }

        /*!
         * \brief Reference to the current solution.
         */
        SpatialFunction &currentSolution()
            { return uCur_; }

        /*!
         * \brief Reference to the right hand side.
         */
        SpatialFunction &rightHandSide()
            { return f_; }

        /*!
         * \brief Reference to solution of the previous time step.
         */
        SpatialFunction &previousSolution()
            { return uPrev_; }

        /*!
         * \brief Reference to solution of the previous time step.
         */
        const SpatialFunction &previousSolution() const
            { return uPrev_; }

        /*!
         * \brief Returns the operator assembler for the global jacobian of
         *        the problem.
         */
        JacobianAssembler &jacobianAssembler()
            { return jacAsm_; }

        /*!
         * \brief Returns the local jacobian which calculates the local
         *        stiffness matrix for an arbitrary cell.
         * 
         * The local stiffness matrices of the cell are used by
         * the jacobian assembler to produce a global linerization of the
         * problem.
         */
        LocalJacobian &localJacobian()
            { return localJacobian_; }

        /*!
         * \brief Same as localJacobian(), included to ease porting.
         */
        LocalJacobian &getLocalJacobian() DUNE_DEPRECATED
            { return localJacobian_; }

        /*!
         * \brief Reference to the grid of the spatial domain.
         */
        const Grid &grid()
            { return problem_.grid(); }

        /*!
         * \brief Try to progress the model to the next timestep.
         */
        template<class NewtonMethod, class NewtonController>
        void update(Scalar &dt, Scalar &nextDt, NewtonMethod &solver, NewtonController &controller)
            {
                asImp_()->updateBegin();

                int numRetries = 0;
                while (true)
                {
                    bool converged = solver.execute(*this->asImp_(),
                                                    controller);
                    nextDt = controller.suggestTimeStepSize(dt);
                    if (converged)
                        break;
                    
                    ++numRetries;
                    if (numRetries > 10)
                        DUNE_THROW(Dune::MathError,
                                   "Newton solver didn't converge after 10 timestep divisions. dt=" << dt);

                    problem_.setTimeStepSize(nextDt);
                    dt = nextDt;
                    
                    asImp_()->updateFailedTry();
                    
                    std::cout << boost::format("Newton didn't converge. Retrying with timestep of %f\n")%dt;
                }
                
                asImp_()->updateSuccessful();
            }

        
        /*!
         * \brief Called by the update() method before it tries to
         *        apply the newton method. This is primary a hook
         *        which the actual model can overload.
         */
        void updateBegin()
            {
                applyDirichletBoundaries_(uCur_);
            }

        
        /*!
         * \brief Called by the update() method if it was
         *        successful. This is primary a hook which the actual
         *        model can overload.
         */
        void updateSuccessful() 
            {
                // make the current solution the previous one.
                *uPrev_ = *uCur_;
            };

        /*!
         * \brief Called by the update() method if a try was
         *         unsuccessful. This is primary a hook which the
         *         actual model can overload.
         */
        void updateFailedTry() 
            {
                // Reset the current solution to the one of the
                // previous time step so that we can start the next
                // update at a physically meaningful solution.
                *uCur_ = *uPrev_;
            };

        /*!
         * \brief Calculate the global residual.
         * 
         * The global deflection of the mass balance from zero.
         */
        void evalGlobalResidual(SpatialFunction &globResidual)
            {
                (*globResidual) = Scalar(0.0);

                // iterate through leaf grid
                CellIterator it     = problem_.grid().template leafbegin<0>();
                CellIterator eendit = problem_.grid().template leafend<0>();
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

                    localJacobian_.setCurrentCell(cell);
                    localJacobian_.restrictToCell(localU, currentSolution());
                    localJacobian_.restrictToCell(localOldU, previousSolution());
                    localJacobian_.setParams(cell, localU, localOldU);
                    localJacobian_.evalLocalResidual(localResidual);

                    // loop over the cell's vertices, map them to the
                    // corresponding grid's node ids and add the
                    // cell's local residual at a node the global
                    // residual at this node.
                    int n = cell.template count<GridDim>();
                    for(int localId=0; localId < n; localId++)
                    {
                        int globalId = problem_.nodeIndex(cell, localId);
                        (*globResidual)[globalId] += localResidual[localId];
                    }
                }
            }


    protected:
        void applyInitialSolution_(SpatialFunction &u)
            {
                // iterate through leaf grid an evaluate c0 at cell center
                CellIterator it     = problem_.grid().template leafbegin<0>();
                CellIterator eendit = problem_.grid().template leafend<0>();
                for (; it != eendit; ++it)
                {
                    // loop over all shape functions of the current cell
                    const Cell& cell = *it;
                    int numNodes = cell.template count<GridDim>();
                    for (int localNodeIdx = 0; localNodeIdx < numNodes; localNodeIdx++) {
                        // get node position in reference coodinates
                        const LocalCoord &local =
                            DomainTraits::referenceElement(it->type()).position(localNodeIdx, GridDim);
                        // get global coordinate of node 
                        const WorldCoord &global = it->geometry()[localNodeIdx];

                        int globalId = this->problem_.nodeIndex(*it, localNodeIdx);
                        
                        // use the problem for actually doing the
                        // dirty work of nailing down the initial
                        // solution.
                        this->problem_.initial((*u)[globalId],
                                               cell,
                                               global,
                                               local);
                    }
                }
            };


        void applyDirichletBoundaries_(SpatialFunction &u)
            {
                // set Dirichlet boundary conditions of the grid's
                // outer boundaries

                UnknownsVector dirichletVal;
                CellIterator cellIt     = problem_.grid().template leafbegin<0>();
                CellIterator cellEndIt  = problem_.grid().template leafend<0>();
                for (; cellIt != cellEndIt; ++cellIt)
                {
                    if (!cellIt->hasBoundaryIntersections())
                        continue;
                    
                    // get the current cell and its set of shape
                    // functions
                    const Cell& cell = *cellIt;
                    Dune::GeometryType geoType = cell.geometry().type();

                    // locally evaluate the cell's boundary condition types
                    localJacobian_.assembleBoundaryCondition(cell);

                    // loop over all the cell's nodes
                    int n = cellIt->template count<GridDim>();
                    for (int i = 0; i < n; ++i) {
                        // translate local node id to a global one
                        int globalId = problem_.nodeIndex(cell, i);

                        // apply dirichlet boundaries but make sure
                        // not to interfere with non-dirichlet
                        // boundaries...
                        bool dirichletEvaluated = false;
                        for (int bcIdx = 0; bcIdx < PrimaryVariables; ++bcIdx) {
                            if (localJacobian_.bc(i)[bcIdx] == BoundaryConditions::dirichlet) {
                                if (!dirichletEvaluated) {
                                    dirichletEvaluated = true;
                                    
                                    // actually evaluate the boundary
                                    // condition for the current cell+node
                                    // combo. 
                                    //
                                    // TODO: better parameters: cell,
                                    //       FVElementGeometry,
                                    //       bfIndex
                                    problem_.dirichlet(dirichletVal,
                                                       cell,
                                                       i,
                                                       globalId);                            
                                }
                                (*u)[globalId][bcIdx] = dirichletVal[bcIdx];
                            }
                        }
                        
                    }
                    
                }
            };

        Implementation *asImp_() 
            { return static_cast<Implementation*>(this); } 
        const Implementation *asImp_() const
            { return static_cast<const Implementation*>(this); } 

        // the problem we want to solve. defines the constitutive
        // relations, material laws, etc.
        Problem     &problem_;

        // the solution we are looking for

        // cur is the current solution, prev the solution of the
        // previous time step
        SpatialFunction uCur_;
        SpatialFunction uPrev_;

        // the right hand side
        SpatialFunction  f_;
        // Linearizes the problem at the current time step using the
        // local jacobian
        JacobianAssembler jacAsm_;
        // calculates the local jacobian matrix for a given cell
        LocalJacobian    &localJacobian_;
    };
}

#endif
