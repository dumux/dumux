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
    typedef typename DomainTraits::Vertex                      Vertex;
    typedef typename DomainTraits::ReferenceElement            ReferenceElement;
    typedef typename BoxTraits::FVElementGeometry              FVElementGeometry;
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

        wasRestarted_ = false;

        // check grid partitioning
        assert((prob.grid().comm().size() == 1) ||
               (prob.grid().overlapSize(0) > 0) ||
               (prob.grid().ghostSize(0) > 0));
    }

    /*!
     * \brief Apply the initial conditions to the model.
     */
    void initial()
    {
        // initialize the static vert data of the box jacobian
        this->localJacobian().setCurSolution(&uCur_);
        this->localJacobian().setOldSolution(&uPrev_);


        if (!wasRestarted_)
        {
            this->localJacobian().initStaticData();
            applyInitialSolution_(uCur_);
        }

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
     * \brief A reference to the problem on which the model is applied.
     */
    const Problem &problem() const
    { return problem_; }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem()
    { return problem_; }

    /*!
     * \brief Reference to the grid of the spatial domain.
     */
    const Grid &grid() const
    { return problem_.grid(); }

    /*!
     * \brief Try to progress the model to the next timestep.
     */
    template<class NewtonMethod, class NewtonController>
    void update(Scalar &dt, Scalar &nextDt, NewtonMethod &solver, NewtonController &controller)
    {
        asImp_().updateBegin();

        int numRetries = 0;
        while (true)
        {
            bool converged = solver.execute(this->asImp_(),
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

            asImp_().updateFailedTry();

            std::cout << boost::format("Newton didn't converge for rank %d. Retrying with timestep of %f\n")
                %grid().comm().rank()%dt;
        }

        asImp_().updateSuccessful();
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
        applyDirichletBoundaries_(uCur_);
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
        ElementIterator        it     = problem_.elementBegin();
        const ElementIterator &eendit = problem_.elementEnd();
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
                int globalId = problem_.vertexIdx(element, localId);
                (*globResidual)[globalId] += localResidual[localId];
            }
        }
    }

    /*!
     * \brief Serializes the current state of the model.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    { res.template serializeEntities<dim>(asImp_(), this->grid()); }

    /*!
     * \brief Deserializes the state of the model.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.template deserializeEntities<dim>(asImp_(), this->grid());
        wasRestarted_ = true;
    }

    /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     */
    void serializeEntity(std::ostream &outstream,
                         const Vertex &vert)
    {
        int vertIdx = problem_.vertexIdx(vert);

        // write phase state
        if (!outstream.good()) {
            DUNE_THROW(IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << (*currentSolution())[vertIdx][eqIdx] << " ";
        }
    };

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     */
    void deserializeEntity(std::istream &instream,
                           const Vertex &vert)
    {
        int vertIdx = problem_.vertexIdx(vert);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(IOError,
                           "Could not deserialize vertex "
                           << vertIdx);
            instream >> (*currentSolution())[vertIdx][eqIdx];
        }
    };

protected:
    void applyInitialSolution_(SpatialFunction &u)
    {
        // first set the whole domain to zero. This is
        // necessary in order to also get a meaningful value
        // for ghost nodes (if we are running in parallel)
        if (problem_.grid().comm().size() > 1) {
            (*u) = Scalar(0.0);
        }

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        //
        // TODO: the initial condition needs to be unique for
        // each vertex. we should think about the API...
        ElementIterator it            = problem_.elementBegin();
        const ElementIterator &eendit = problem_.elementEnd();
        for (; it != eendit; ++it)
        {
            // deal with the current element
            applyInitialSolutionElement_(u, *it);
        }
    };

    // apply the initial solition for a single element
    void applyInitialSolutionElement_(SpatialFunction &u,
                                      const Element &element)
    {
        // HACK: set the current element for the local
        // solution in order to get an updated FVElementGeometry
        localJacobian_.setCurrentElement(element);

        // loop over all element vertices, i.e. sub control volumes
        int numScv = element.template count<dim>();
        for (int scvIdx = 0;
             scvIdx < numScv;
             scvIdx++)
        {
            int globalIdx = this->problem_.vertexIdx(element,
                                                     scvIdx);

            const FVElementGeometry &fvElemGeom
                = localJacobian_.curFvElementGeometry();
            // use the problem for actually doing the
            // dirty work of nailing down the initial
            // solution.
            this->problem_.initial((*u)[globalIdx],
                                   element,
                                   fvElemGeom,
                                   scvIdx);
        }
    }


    // apply dirichlet boundaries for the whole grid
    void applyDirichletBoundaries_(SpatialFunction &u)
    {
        // set Dirichlet boundary conditions of the grid's
        // outer boundaries
        ElementIterator elementIt           = problem_.elementBegin();
        const ElementIterator &elementEndIt = problem_.elementEnd();
        for (; elementIt != elementEndIt; ++elementIt)
        {
            // ignore elements which are not on the boundary of
            // the domain
            if (!elementIt->hasBoundaryIntersections())
                continue;

            // evaluate the element's boundary locally
            localJacobian_.updateBoundaryTypes(*elementIt);

            // apply dirichlet boundary for the current element
            applyDirichletElement_(u, *elementIt);
        }
    };

    // apply dirichlet boundaries for a single element
    void applyDirichletElement_(SpatialFunction &u,
                                const Element &element)
    {
        Dune::GeometryType      geoType = element.geometry().type();
        const ReferenceElement &refElem = DomainTraits::referenceElement(geoType);

        // loop over all the element's surface patches
        IntersectionIterator isIt = element.ileafbegin();
        const IntersectionIterator &isEndIt = element.ileafend();
        for (; isIt != isEndIt; ++isIt) {
            // if the current intersection is not on the boundary,
            // we ignore it
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the
            // current face
            int faceIdx = isIt->numberInInside();
            int numVerticesOfFace = refElem.size(faceIdx, 1, dim);
            for (int vertInFace = 0;
                 vertInFace < numVerticesOfFace;
                 vertInFace++)
            {
                // apply dirichlet boundaries for the current
                // sub-control volume face
                applyDirichletSCVF_(u, element, refElem, isIt, vertInFace);
            }
        }
    }

    // apply dirichlet boundaries for a single boundary
    // sub-control volume face of a finite volume cell.
    void applyDirichletSCVF_(SpatialFunction &u,
                             const Element &element,
                             const ReferenceElement &refElem,
                             const IntersectionIterator &isIt,
                             int faceVertIdx)
    {
        // apply dirichlet boundaries but make sure
        // not to interfere with non-dirichlet
        // boundaries...
        const FVElementGeometry &fvElemGeom =
            localJacobian_.curFvElementGeometry();

        int elemVertIdx = refElem.subEntity(isIt->numberInInside(),
                                            1,
                                            faceVertIdx,
                                            dim);
        int boundaryFaceIdx = fvElemGeom.boundaryFaceIndex(isIt->numberInInside(),
                                                           faceVertIdx);
        int globalVertexIdx = problem_.vertexIdx(element,
                                                 elemVertIdx);


        SolutionVector dirichletVal;
        bool dirichletEvaluated = false;
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            // ignore non-dirichlet boundary conditions
            if (localJacobian_.bc(elemVertIdx)[eqIdx] != BoundaryConditions::dirichlet)
                continue;

            // make sure to evaluate the dirichlet boundary
            // conditions exactly once (and only if the boundary
            // type is actually dirichlet).
            if (!dirichletEvaluated)
            {
                dirichletEvaluated = true;
                problem_.dirichlet(dirichletVal,
                                   element,
                                   fvElemGeom,
                                   isIt,
                                   elemVertIdx,
                                   boundaryFaceIdx);
            }

            // copy the dirichlet value for the current equation
            // to the global index.
            //
            // TODO: we should use the sum weighted by the
            //       sub-control volume face's area instead of
            //       just overwriting the previous values...
            (*u)[globalVertexIdx][eqIdx] = dirichletVal[eqIdx];
        }
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

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

    bool wasRestarted_;
};
}

#endif
