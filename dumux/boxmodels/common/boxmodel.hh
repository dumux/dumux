// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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

#include <dumux/common/valgrind.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <boost/format.hpp>

#include "boxproperties.hh"

#include "boxelementvolumevariables.hh"
#include "boxlocaljacobian.hh"
#include "boxlocalresidual.hh"

#include "pdelabboxlocaloperator.hh"

namespace Dumux
{

/*!
 * \defgroup BoxModel   Box-Scheme
 */
/*!
 * \ingroup BoxModel
 * \defgroup BoxProblems Box-Problems
 */
/*!
 * \ingroup BoxModel
 * \defgroup BoxModels   Box-Models
 */


/*!
 * \ingroup BoxModel
 *
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class BoxModel
{
    typedef BoxModel<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DofMapper)) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianAssembler)) JacobianAssembler;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        enableJacobianRecycling  = GET_PROP_VALUE(TypeTag, PTAG(EnableJacobianRecycling)),
        dim = GridView::dimension
    };

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;
    typedef typename RefElemProp::ReferenceElement ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

public:
    /*!
     * \brief The constructor.
     */
    BoxModel()
    { }

    ~BoxModel()
    { delete jacAsm_;  }

    /*!
     * \brief Apply the initial conditions to the model.
     */
    void init(Problem &prob)
    {
        problemPtr_ = &prob;

        int nDofs = asImp_().numDofs();
        uCur_.resize(nDofs);
        uPrev_.resize(nDofs);
        boxVolume_.resize(nDofs);

        localJacobian_.init(problem_());
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(problem_());

        applyInitialSolution_();

        // also set the solution of the "previous" time step to the
        // initial solution.
        uPrev_ = uCur_;
    }

    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     */
    Scalar globalResidual(SolutionVector &dest, const SolutionVector &u)
    {
        SolutionVector tmp(curSol());
        curSol() = u;
        Scalar res = globalResidual(dest);
        curSol() = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     */
    Scalar globalResidual(SolutionVector &dest)
    {
        dest = 0;

        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            localResidual().eval(*elemIt);

            for (int i = 0; i < elemIt->template count<dim>(); ++i) {
                int globalI = vertexMapper().map(*elemIt, i, dim);
                dest[globalI] += localResidual().residual(i);
            }
        };

        Scalar result = dest.two_norm();
        /*
        Scalar result = 0;
        for (int i = 0; i < (*tmp).size(); ++i) {
            for (int j = 0; j < numEq; ++j)
                result += std::abs((*tmp)[i][j]);
        }
        */
        return result;
    }

    /*!
     * \brief Compute the total storage of all conservation quantities.
     */
    void globalStorage(PrimaryVariables &dest)
    {
        dest = 0;
        
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            localResidual().evalStorage(*elemIt);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                dest += localResidual().residual(i);
        };

        gridView_().comm().sum(dest);
    }

    /*!
     * \brief Returns the volume of a given control volume.
     */
    Scalar boxVolume(int globalIdx) const
    { return boxVolume_[globalIdx][0]; }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    const SolutionVector &curSol() const
    { return uCur_; }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    SolutionVector &curSol()
    { return uCur_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    const SolutionVector &prevSol() const
    { return uPrev_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    SolutionVector &prevSol()
    { return uPrev_; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

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
     * \copydoc localJacobian()
     */
    const LocalJacobian &localJacobian() const
    { return localJacobian_; }

    /*!
     * \brief Returns the local residual function.
     */
    LocalResidual &localResidual()
    { return localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual()
     */
    const LocalResidual &localResidual() const
    { return localJacobian().localResidual(); }

    /*!
     * \brief Try to progress the model to the next timestep.
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i)
            Valgrind::CheckDefined(curSol()[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.execute(controller);
        if (converged)
            asImp_().updateSuccessful();
        else
            asImp_().updateFailed();

#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i) {
            Valgrind::CheckDefined(curSol()[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }


    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primary a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    { }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    {
        // the jacobian matrix of the last iteration of the current
        // time step is the same as the jacobian of the first
        // iteration of the next time step
        if (enableJacobianRecycling)
            jacAsm_->setMatrixReuseable(true);

        // make the current solution the previous one.
        uPrev_ = uCur_;
    };

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        uCur_ = uPrev_;
    };

    /*!
     * \brief Serializes the current state of the model.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    { res.template serializeEntities<dim>(asImp_(), this->gridView_()); }

    /*!
     * \brief Deserializes the state of the model.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.template deserializeEntities<dim>(asImp_(), this->gridView_());
        prevSol() = curSol();
    }

    /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     */
    void serializeEntity(std::ostream &outstream,
                         const Vertex &vert)
    {
        int vertIdx = dofMapper().map(vert);

        // write phase state
        if (!outstream.good()) {
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << curSol()[vertIdx][eqIdx] << " ";
        }
    };

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     */
    void deserializeEntity(std::istream &instream,
                           const Vertex &vert)
    {
        int vertIdx = dofMapper().map(vert);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << vertIdx);
            instream >> curSol()[vertIdx][eqIdx];
        }
    };

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    { return gridView_().size(dim); }

    /*!
     * \brief Mapper for the entities where degrees of freedoms are
     *        defined to indices.
     *
     * This usually means a mapper for vertices.
     */
    const DofMapper &dofMapper() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); };

    void resetJacobianAssembler ()
    {
        delete jacAsm_;
        jacAsm_ = new JacobianAssembler;
        jacAsm_->init(problem_());
    }

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK multi writer.
     *
     * \param writer  The VTK multi writer where the fields should be added.
     * \param update  The delte of the solution function before and after the Newton update
     */
    template <class MultiWriter>
    void addConvergenceVtkFields(MultiWriter &writer,
                                 const SolutionVector &u,
                                 const SolutionVector &deltaU)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        SolutionVector globalResid(asImp_(), 0.0);
        asImp_().globalResidual(globalResid, u);

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        //unsigned numElements = this->gridView_().size(0);

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        for (int i = 0; i < numEq; ++i) {
            x[i] = writer.template createField<Scalar, 1>(numVertices);
            delta[i] = writer.template createField<Scalar, 1>(numVertices);
            def[i] = writer.template createField<Scalar, 1>(numVertices);
        }

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++ vIt)
        {
            int globalIdx = vertexMapper().map(*vIt);
            for (int i = 0; i < numEq; ++i) {
                (*x[i])[globalIdx] = u[globalIdx][i];
                (*delta[i])[globalIdx] = - deltaU[globalIdx][i];
                (*def[i])[globalIdx] = globalResid[globalIdx][i];
            }
        }

        for (int i = 0; i < numEq; ++i) {
            writer.addVertexData(x[i], (boost::format("x_%i")%i).str().c_str());
            writer.addVertexData(delta[i], (boost::format("delta_%i")%i).str().c_str());
            writer.addVertexData(def[i], (boost::format("defect_%i")%i).str().c_str());
        }

        asImp_().addOutputVtkFields(u, writer);
    }

    /*!
     * \brief Add the quantities of a time step which ought to be
     *        written to disk.
     *
     * This should be overwritten by the acutal model if any secondary
     * variables should be written out. Read: This should _always_ be
     * overwritten by well behaved models!
     *
     * \param writer  The VTK multi writer where the fields should be added.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);

        // global defect of the two auxiliary equations
        ScalarField* x[numEq];
        for (int i = 0; i < numEq; ++i) {
            x[i] = writer.template createField<Scalar, 1>(numVertices);
        }

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++ vIt)
        {
            int globalIdx = vertexMapper().map(*vIt);
            for (int i = 0; i < numEq; ++i) {
                (*x[i])[globalIdx] = sol[globalIdx][i];
            }
        }

        for (int i = 0; i < numEq; ++i)
            writer.addVertexData(x[i], (boost::format("primaryVar%i")%i).str().c_str());
    }


protected:
    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    LocalResidual &localResidual_()
    { return localJacobian_.localResidual(); }

    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        uCur_ = Scalar(0.0);
        boxVolume_ = Scalar(0.0);

        FVElementGeometry fvElemGeom;

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        //
        // TODO: the initial condition needs to be unique for
        // each vertex. we should think about the API...
        ElementIterator it = gridView_().template begin<0>();
        const ElementIterator &eendit = gridView_().template end<0>();
        for (; it != eendit; ++it) {
            // deal with the current element
            fvElemGeom.update(gridView_(), *it);

            // loop over all element vertices, i.e. sub control volumes
            for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; scvIdx++)
            {
                // map the local vertex index to the global one
                int globalIdx = vertexMapper().map(*it,
                                                   scvIdx,
                                                   dim);

                // let the problem do the dirty work of nailing down
                // the initial solution.
                PrimaryVariables initVal;
                Valgrind::SetUndefined(initVal);
                problem_().initial(initVal,
                                   *it,
                                   fvElemGeom,
                                   scvIdx);
                Valgrind::CheckDefined(initVal);

                // add up the initial values of all sub-control
                // volumes. If the initial values disagree for
                // different sub control volumes, the initial value
                // will be the arithmetic mean.
                initVal *= fvElemGeom.subContVol[scvIdx].volume;
                boxVolume_[globalIdx] += fvElemGeom.subContVol[scvIdx].volume;
                uCur_[globalIdx] += initVal;
                Valgrind::CheckDefined(uCur_[globalIdx]);
            }
        }
        // divide all primary variables by the volume of their boxes
        int n = gridView_().size(dim);
        for (int i = 0; i < n; ++i) {
            uCur_[i] /= boxVolume(i);
        }
    }

    bool verbose_() const
    { return gridView_().comm().rank() == 0; };

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    JacobianAssembler *jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;
};
}

#endif
