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
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Caculates the Jacobian of the local residual for box models
 */
#ifndef DUMUX_BOX_LOCAL_JACOBIAN_HH
#define DUMUX_BOX_LOCAL_JACOBIAN_HH

#include <dune/istl/matrix.hh>

#include "boxelementboundarytypes.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Caculates the Jacobian of the local residual for box models
 *
 * The default behaviour is to use numeric differentiation, i.e.
 * forward or backward differences (2nd order), or central
 * differences (3rd order). The method used is determined by the
 * "NumericDifferenceMethod" property:
 *
 * - if the value of this property is smaller than 0, backward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
 *   \f]
 *
 * - if the value of this property is 0, central
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
 *   \f]
 *
 * - if the value of this property is larger than 0, forward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
 *   \f]
 *
 * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
 * is the value of a sub-control volume's primary variable at the
 * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianAssembler)) JacobianAssembler;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        Red = JacobianAssembler::Red,
        Yellow = JacobianAssembler::Yellow,
        Green = JacobianAssembler::Green,

        numDiffMethod = GET_PROP_VALUE(TypeTag,
                                       PTAG(NumericDifferenceMethod))
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Element::EntityPointer ElementPointer;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    // copying a local jacobian is not a good idea
    BoxLocalJacobian(const BoxLocalJacobian &);

public:
    BoxLocalJacobian()
    { Valgrind::SetUndefined(problemPtr_); }


    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param prob The problem which we want to simulate.
     */
    void init(Problem &prob)
    {
        problemPtr_ = &prob;
        localResidual_.init(prob);
        // assume quadrilinears as elements with most vertices
        A_.setSize(2<<dim, 2<<dim);
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    void assemble(const Element &element)
    {
        // set the current grid element and update the element's
        // finite volume geometry
        elemPtr_ = &element;
        fvElemGeom_.update(gridView_(), element);
        reset_();

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
        // calculate the local residual
        localResidual().eval(elem_(),
                             fvElemGeom_,
                             prevVolVars_,
                             curVolVars_,
                             bcTypes_);
        residual_ = localResidual().residual();

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
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    /*!
     * \brief Returns the Jacobian of the equations at vertex i to the
     *        primary variables at vertex j.
     *
     * \param i The local vertex (or sub-contol volume) index on which
     *          the equations are defined
     * \param j The local vertex (or sub-contol volume) index which holds
     *          primary variables
     */
    const MatrixBlock &mat(int i, int j) const
    { return A_[i][j]; }

    /*!
     * \brief Returns the residual of the equations at vertex i.
     *
     * \param i The local vertex (or sub-contol volume) index on which
     *          the equations are defined
     */
    const PrimaryVariables &residual(int i) const
    { return residual_[i]; }

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
     * \brief Returns a reference to the jacobian assembler.
     */
    const JacobianAssembler &jacAsm_() const
    { return model_().jacobianAssembler(); }

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Reset the local jacobian matrix to 0
     */
    void reset_()
    {
        int n = elem_().template count<dim>();
        for (int i = 0; i < n; ++ i) {
            for (int j = 0; j < n; ++ j) {
                A_[i][j] = 0.0;
            }
        }
    }

    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than numerical differentiation is available.
     *
     * The default implementation of this method uses numeric
     * differentiation, i.e. forward or backward differences (2nd
     * order), or central differences (3rd order). The method used is
     * determined by the "NumericDifferenceMethod" property:
     *
     * - if the value of this property is smaller than 0, backward
     *   differences are used, i.e.:
     *   \f[
         \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
     *   \f]
     *
     * - if the value of this property is 0, central
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
     *   \f]
     *
     * - if the value of this property is larger than 0, forward
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
     *   \f]
     *
     * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
     * is the value of a sub-control volume's primary variable at the
     * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
     *
     * \param dest The vector storing the partial derivatives of all
     *              equations
     * \param scvIdx The sub-control volume index of the current
     *               finite element for which the partial derivative
     *               ought to be calculated
     * \param pvIdx The index of the primary variable at the scvIdx'
     *              sub-control volume of the current finite element
     *              for which the partial derivative ought to be
     *              calculated
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
        Scalar delta = 0;

        if (numDiffMethod >= 0) {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

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
        }
        else {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            dest = residual_;
        }


        if (numDiffMethod <= 0) {
            // we are not using forward differences, i.e. we don't
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

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
            dest -= localResidual().residual();
        }
        else {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            dest -= residual_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        dest /= delta;

        // restore the orignal state of the element's volume variables
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
     * \brief Updates the current local Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at vertex 'scvIdx' .
     */
    void updateLocalJacobian_(int scvIdx,
                              int pvIdx,
                              const ElementSolutionVector &deriv)
    {
        for (int i = 0; i < fvElemGeom_.numVertices; i++)
        {
            if (jacAsm_().vertexColor(elem_(), i) == Green) {
                // Green vertices are not to be changed!
                continue;
            }

            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[i][scvIdx][eqIdx][pvIdx] is the rate of change of
                // the residual of equation 'eqIdx' at vertex 'i'
                // depending on the primary variable 'pvIdx' at vertex
                // 'scvIdx'.
                this->A_[i][scvIdx][eqIdx][pvIdx] = deriv[i][eqIdx];
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
    ElementSolutionVector residual_;
};
}

#endif
