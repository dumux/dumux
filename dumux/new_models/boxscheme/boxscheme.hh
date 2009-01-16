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
        typedef typename DomainTraits::Element                     Element;
        typedef typename DomainTraits::ReferenceElement            ReferenceElement;
        typedef typename DomainTraits::ElementIterator             ElementIterator;
        typedef typename DomainTraits::IntersectionIterator        IntersectionIterator;
        typedef typename DomainTraits::CoordScalar                 CoordScalar;
        typedef typename DomainTraits::GlobalPosition              GlobalPosition;
        typedef typename DomainTraits::LocalPosition               LocalPosition;

        typedef typename BoxTraits::JacobianAssembler          JacobianAssembler;
        typedef typename BoxTraits::SpatialFunction            SpatialFunction;
        typedef typename SpatialFunction::RepresentationType   BoxFnRep;
        typedef typename BoxTraits::LocalFunction              LocalFunction;

        typedef typename BoxTraits::ShapeFunctionSetContainer  ShapeFunctionSetContainer;
        typedef typename ShapeFunctionSetContainer::value_type ShapeFunctionSet;

        typedef typename BoxTraits::SolutionVector      SolutionVector;
        typedef typename BoxTraits::BoundaryTypeVector  BoundaryTypeVector;

        typedef LocalJacobianT                          LocalJacobian;

        // some constants
        enum {
            numEq = BoxTraits::numEq,

            dim      = DomainTraits::dim,
            dimWorld = DomainTraits::dimWorld
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

                // check data partitioning
                assert((prob.grid().overlapSize(0) > 0) || (prob.grid().ghostSize(0) > 0));
            }

        /*!
         * \brief Apply the initial conditions to the model.
         */
        void initial()
            {
                // initialize the static vert data of the box jacobian
                this->localJacobian().setCurSolution(&uCur_);
                this->localJacobian().setOldSolution(&uPrev_);

                this->localJacobian().initStaticData();

                applyInitialSolution_(uCur_);
                applyDirichletBoundaries_(uCur_);

                *uPrev_ = *uCur_;

                // update the static vert data with the initial solution
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
         *        stiffness matrix for an arbitrary element.
         *
         * The local stiffness matrices of the element are used by
         * the jacobian assembler to produce a global linerization of the
         * problem.
         */
        LocalJacobian &localJacobian()
            { return localJacobian_; }

        /*!
         * \brief Reference to the grid of the spatial domain.
         */
        const Grid &grid() const
            { return problem_.grid(); }

        /*!
         * \brief Reference to the grid of the spatial domain.
         */
/*        Grid &grid()
            { return problem_.grid(); }
*/

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
                    if (converged) {
                        std::cout << boost::format("Newton solver converged for rank %d\n")
                            %grid().comm().rank();
                        break;
                    }

                    ++numRetries;
                    if (numRetries > 10)
                        DUNE_THROW(Dune::MathError,
                                   "Newton solver didn't converge after 10 timestep divisions. dt=" << dt);

                    problem_.setTimeStepSize(nextDt);
                    dt = nextDt;

                    asImp_()->updateFailedTry();

                    std::cout << boost::format("Newton didn't converge for rank %d. Retrying with timestep of %f\n")
                        %grid().comm().rank()%dt;
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
                ElementIterator it     = problem_.grid().template leafbegin<0>();
                ElementIterator eendit = problem_.grid().template leafend<0>();
                for (; it != eendit; ++it)
                {
                    // tell the local jacobian which element it should
                    // consider and evaluate the local residual for the
                    // element. in order to do this we first have to
                    // evaluate the element's local solutions for the
                    // current and the last timestep.
                    const Element& element = *it;
                    const int numVertices = element.template count<dim>();
                    LocalFunction localResidual(numVertices);
                    LocalFunction localU(numVertices);
                    LocalFunction localOldU(numVertices);

                    localJacobian_.setCurrentElement(element);
                    localJacobian_.restrictToElement(localU, currentSolution());
                    localJacobian_.restrictToElement(localOldU, previousSolution());
                    localJacobian_.setParams(element, localU, localOldU);
                    localJacobian_.evalLocalResidual(localResidual);

                    // loop over the element's vertices, map them to the
                    // corresponding grid's vert ids and add the
                    // element's local residual at a vert the global
                    // residual at this vert.
                    int n = element.template count<dim>();
                    for(int localId=0; localId < n; localId++)
                    {
                        int globalId = problem_.vertIdx(element, localId);
                        (*globResidual)[globalId] += localResidual[localId];
                    }
                }
            }


    protected:
        void applyInitialSolution_(SpatialFunction &u)
            {
                // iterate through leaf grid an evaluate c0 at element center
                ElementIterator it     = problem_.grid().template leafbegin<0>();
                ElementIterator eendit = problem_.grid().template leafend<0>();
                for (; it != eendit; ++it)
                {
                    // loop over all shape functions of the current element
                    const Element& element = *it;
                    int numVertices = element.template count<dim>();
                    for (int localVertexIdx = 0; localVertexIdx < numVertices; localVertexIdx++) {
                        // get vert position in reference coodinates
                        const LocalPosition &local =
                            DomainTraits::referenceElement(it->type()).position(localVertexIdx, dim);
                        // get global coordinate of vert
                        const GlobalPosition &global = it->geometry().corner(localVertexIdx);

                        int globalId = this->problem_.vertIdx(*it, localVertexIdx);

                        // use the problem for actually doing the
                        // dirty work of nailing down the initial
                        // solution.
                        this->problem_.initial((*u)[globalId],
                                               element,
                                               global,
                                               local);
                    }
                }
            };


        void applyDirichletBoundaries_(SpatialFunction &u)
            {
                // set Dirichlet boundary conditions of the grid's
                // outer boundaries
                SolutionVector dirichletVal(0);
                ElementIterator elementIt     = problem_.grid().template leafbegin<0>();
                ElementIterator elementEndIt  = problem_.grid().template leafend<0>();
                for (; elementIt != elementEndIt; ++elementIt)
                {
                    if (!elementIt->hasBoundaryIntersections())
                        continue;

                    // get the current element and its set of shape
                    // functions
                    const Element& element = *elementIt;
                    Dune::GeometryType geoType = element.geometry().type();

                    // locally evaluate the element's boundary condition types
                    localJacobian_.assembleBoundaryCondition(element);

                    // loop over all the element's verts
                    int n = elementIt->template count<dim>();
                    for (int i = 0; i < n; ++i) {
                        // translate local vert id to a global one
                        int globalId = problem_.vertIdx(element, i);

                        // apply dirichlet boundaries but make sure
                        // not to interfere with non-dirichlet
                        // boundaries...
                        bool dirichletEvaluated = false;
                        for (int bcIdx = 0; bcIdx < numEq; ++bcIdx) {
                            if (localJacobian_.bc(i)[bcIdx] == BoundaryConditions::dirichlet) {
                                if (!dirichletEvaluated) {
                                    dirichletEvaluated = true;

                                    // actually evaluate the boundary
                                    // condition for the current element+vert
                                    // combo.
                                    //
                                    // TODO: better parameters: element,
                                    //       FVElementGeometry,
                                    //       bfIdx
                                    problem_.dirichlet(dirichletVal,
                                                       element,
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
        // calculates the local jacobian matrix for a given element
        LocalJacobian    &localJacobian_;
    };
}

#endif
