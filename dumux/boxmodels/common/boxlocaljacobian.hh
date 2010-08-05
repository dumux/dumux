// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
/*!
 * \file
 * \brief Caculates the jacobian of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_LOCAL_JACOBIAN_HH
#define DUMUX_BOX_LOCAL_JACOBIAN_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/valgrind.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <boost/format.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/matrix.hh>

#include "boxelementvolumevariables.hh"
#include "boxfvelementgeometry.hh"
#include "boxlocalresidual.hh"

#include "boxproperties.hh"



namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Element-wise caculation of the jacobian matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalJacobian
{
private:
    typedef BoxLocalJacobian<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Element::EntityPointer ElementPointer;

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;
    typedef typename RefElemProp::ReferenceElement ReferenceElement;

    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

public:
    BoxLocalJacobian()
    {
        Valgrind::SetUndefined(problemPtr_);
    }


    void init(Problem &prob)
    {
        problemPtr_ = &prob;
        localResidual_.init(prob);
        // assume quadrilinears as elements with most vertices
        A_.setSize(2<<dim, 2<<dim);
    }

    /*!
     * \brief Assemble the linear system of equations for the
     *        verts of a element, given a local solution 'localU'.
     */
    void assemble(const Element &element)
    {
        // set the current grid element and update the element's
        // finite volume geometry
        elemPtr_ = &element;
        fvElemGeom_.update(gridView_(), elem_());
        bcTypes_.update(problem_(), elem_(), fvElemGeom_);

        // this is pretty much a HACK because the internal state of
        // the problem is not supposed to be changed during the
        // evaluation of the residual. (Reasons: It is a violation of
        // abstraction, makes everything more prone to errors and is
        // not thread save.) The real solution are context objects!
        problem_().updateCouplingParams(elem_());


        int numVertices = fvElemGeom_.numVertices;

        // update the secondary variables for the element at the last
        // and the current time levels
        prevVolVars_.update(problem_(),
                            elem_(),
                            fvElemGeom_,
                            true /* isOldSol? */);
        curVolVars_.update(problem_(),
                           elem_(),
                           fvElemGeom_,
                           false /* isOldSol? */);

        // calculate the local jacobian matrix
        ElementSolutionVector partialDeriv(numVertices);
        for (int j = 0; j < numVertices; j++) {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().evalPartialDerivative_(partialDeriv,
                                                j,
                                                pvIdx);

                // update the local stiffness matrix with the current partial
                // derivatives
                updateLocalJacobian_(j,
                                     pvIdx,
                                     partialDeriv);
            }
        }
    }

    /*!
     * \brief Returns a reference to the local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns a reference to the local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    /*!
     * \brief Returns the jacobian of the equations at vertex i to the
     *        primary variables at vertex j.
     */
    const MatrixBlock &mat(int i, int j) const
    { return A_[i][j]; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    {
        Valgrind::CheckDefined(problemPtr_);
        return *problemPtr_;
    };

    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the element.
     */
    const Element &elem_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    };

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); };

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than central differences ought to be used.
     */
    void evalPartialDerivative_(ElementSolutionVector &dest,
                                int scvIdx,
                                int pvIdx)
    {
        int globalIdx = vertexMapper_().map(elem_(), scvIdx, dim);
        PrimaryVariables priVars(model_().curSol()[globalIdx]);
        VolumeVariables origVolVars(curVolVars_[scvIdx]);

        curVolVars_[scvIdx].setEvalPoint(&origVolVars);
        Scalar eps = asImp_().numericEpsilon_(scvIdx, pvIdx);

        // deflect primary variables
        priVars[pvIdx] += eps;

        // calculate the residual
        curVolVars_[scvIdx].update(priVars,
                                   problem_(),
                                   elem_(),
                                   fvElemGeom_,
                                   scvIdx,
                                   false);
        localResidual().eval(elem_(),
                             fvElemGeom_,
                             prevVolVars_,
                             curVolVars_,
                             bcTypes_);

        // store the residual
        dest = localResidual().residual();

        // deflect the primary variables
        priVars[pvIdx] -= 2*eps;

        // calculate residual again
        curVolVars_[scvIdx].update(priVars,
                                   problem_(),
                                   elem_(),
                                   fvElemGeom_,
                                   scvIdx,
                                   false);
        localResidual().eval(elem_(),
                             fvElemGeom_,
                             prevVolVars_,
                             curVolVars_,
                             bcTypes_);

        // central differences
        dest -= localResidual().residual();
        dest /= 2*eps;

        // restore the orignal state of the element's secondary
        // variables
        curVolVars_[scvIdx] = origVolVars;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < dest.size(); ++i)
            Valgrind::CheckDefined(dest[i]);
#endif
    }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemSol    The current solution on the element
     * \param scvIdx     The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon_(int scvIdx,
                           int pvIdx) const
    {
        Scalar pv = this->curVolVars_[scvIdx].primaryVars()[pvIdx];
        return 1e-9*(std::abs(pv) + 1);
    }

    /*!
     * \brief Updates the current local stiffness matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at vertex j .
     */
    void updateLocalJacobian_(int scvIdx,
                              int pvIdx,
                              const ElementSolutionVector &stiffness)
    {
        for (int i = 0; i < fvElemGeom_.numVertices; i++) {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[i][scvIdx][eqIdx][pvIdx] is the approximate rate
                // of change of the residual of equation 'eqIdx' at
                // vertex 'i' depending on the primary variable
                // 'pvIdx' at vertex 'scvIdx'.
                this->A_[i][scvIdx][eqIdx][pvIdx] = stiffness[i][eqIdx];
                Valgrind::CheckDefined(this->A_[i][scvIdx][eqIdx][pvIdx]);
            }
        }
    }

    const Element *elemPtr_;

    FVElementGeometry fvElemGeom_;
    ElementBoundaryTypes bcTypes_;

    // The problem we would like to solve
    Problem *problemPtr_;

    // secondary variables at the previous and at the current time
    // levels
    ElementVolumeVariables prevVolVars_;
    ElementVolumeVariables curVolVars_;

    LocalResidual localResidual_;
    LocalBlockMatrix A_;
};
}

#endif
