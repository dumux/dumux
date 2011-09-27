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

#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>

#include "boxelementboundarytypes.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Calculates the Jacobian of the local residual for box models
 *
 * The default behavior is to use numeric differentiation, i.e.
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

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVariables)) ElementVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;
    typedef Dune::BlockVector<MatrixBlock> LocalStorageMatrix;

    // copying a local jacobian is not a good idea
    BoxLocalJacobian(const BoxLocalJacobian &);

public:
    BoxLocalJacobian()
    { 
        internalElemVars_ = 0;
    }

    ~BoxLocalJacobian()
    {
        delete internalElemVars_;
    }

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
        modelPtr_ = &prob.model();
        internalElemVars_ = new ElementVariables(prob);
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * This assembles the 'grad f(x^k)' and 'f(x^k)' part of the newton update
     */
    void assemble(const Element &element)
    {
        internalElemVars_->updateAll(element);

        assemble(*internalElemVars_);
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect, given all secondary variables for the element.
     *
     * After calling this method the ElementVariables are in undefined
     * state, so do not use it anymore!
     */
    void assemble(ElementVariables &elemVars)
    {        
        // update the weights of the primary variables using the
        // current element variables
        model_().updatePVWeights(elemVars);
        
        resize_(elemVars);
        reset_(elemVars);
       
        // calculate the local residual
        localResidual_.eval(residual_, residualStorage_, elemVars);

        // save all flux variables calculated using the unmodified
        // primary variables. This automatically makes these flux
        // variables the evaluation point.
        elemVars.saveScvfVars();

        // calculate the local jacobian matrix
        int numScv = elemVars.numScv();
        for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().evalPartialDerivative_(elemVars,
                                                scvIdx,
                                                pvIdx);

                // update the local stiffness matrix with the current
                // partial derivatives
                updateLocalJacobian_(elemVars, scvIdx, pvIdx);
            }
        }

        // restore flux variables. 
        //elemVars.restoreScvfVars(); // not necessary
    }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param scvIdx     The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon(const ElementVariables &elemVars, 
                          int scvIdx,
                          int pvIdx) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        static const Scalar baseEps 
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(),
                                           1.0);
        
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        Scalar pv = elemVars.volVars(scvIdx, /*historyIdx=*/0).primaryVar(pvIdx);
        return baseEps*(std::abs(pv) + 1);
    }

    /*!
     * \brief Return reference to the local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns the local Jacobian matrix of the residual of a sub-control volume.
     *
     * \param domainScvIdx The local index of the sub control volume which contains the independents
     * \param rangeScvIdx The local index of the sub control volume which contains the local residual
     */
    const MatrixBlock &jacobian(int domainScvIdx, int rangeScvIdx) const
    { return jacobian_[domainScvIdx][rangeScvIdx]; }

    /*!
     * \brief Returns the local Jacobian matrix the storage term of a sub-control volume.
     *
     * \param scvIdx The local index of sub control volume
     */
    const MatrixBlock &jacobianStorage(int scvIdx) const
    { return jacobianStorage_[scvIdx]; }

    /*!
     * \brief Returns the local residual of a sub-control volume.
     *
     * \param scvIdx The local index of the sub control volume
     */
    const VectorBlock &residual(int scvIdx) const
    { return residual_[scvIdx]; }

    /*!
     * \brief Returns the local storage term of a sub-control volume.
     *
     * \param scvIdx The local index of the sub control volume
     */
    const VectorBlock &residualStorage(int scvIdx) const
    { return residualStorage_[scvIdx]; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Problem &problem_() const
    { return *problemPtr_; }
    const Model &model_() const
    { return *modelPtr_; }

    /*!
     * \brief Returns the numeric difference method which is applied.
     */
    static int numericDifferenceMethod_() 
    { return GET_PARAM(TypeTag, int, NumericDifferenceMethod); }

    /*!
     * \brief Resize all internal attributes to the size of the
     *        element.
     */
    void resize_(const ElementVariables &elemVars)
    {
        int n = elemVars.numScv();

        jacobian_.setSize(n, n);
        jacobianStorage_.resize(n);
        
        residual_.resize(n);
        residualStorage_.resize(n);
        
        derivResidual_.resize(n);
        derivStorage_.resize(n);
    };

    /*!
     * \brief Reset the all relevant internal attributes to 0
     */
    void reset_(const ElementVariables &elemVars)
    {
        int numScv = elemVars.numScv();
        for (int i = 0; i < numScv; ++ i) {
            residual_[i] = 0.0;
            residualStorage_[i] = 0.0;

            jacobianStorage_[i] = 0.0;
            for (int j = 0; j < numScv; ++ j) {
                jacobian_[i][j] = 0.0;
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
    void evalPartialDerivative_(ElementVariables &elemVars,
                                int scvIdx,
                                int pvIdx)
    {
        // save all quantities which depend on the specified primary
        // variable at the given sub control volume
        elemVars.saveScvVars(scvIdx);

        PrimaryVariables priVars(elemVars.volVars(scvIdx, /*historyIdx=*/0).primaryVars());
        Scalar eps = asImp_().numericEpsilon(elemVars, scvIdx, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_() >= 0) {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // calculate the residual
            elemVars.updateScvVars(priVars, scvIdx, /*historyIdx=*/0);
            elemVars.updateAllScvfVars();
            localResidual_.eval(derivResidual_, derivStorage_, elemVars);
        }
        else {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            derivResidual_ = residual_;
            derivStorage_[scvIdx] = residualStorage_[scvIdx];
        }

        if (numericDifferenceMethod_() <= 0) {
            // we are not using forward differences, i.e. we don't
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // calculate residual again, this time we use the local
            // residual's internal storage.
            elemVars.updateScvVars(priVars, scvIdx, /*historyIdx=*/0);
            elemVars.updateAllScvfVars();
            localResidual_.eval(elemVars);
            
            derivResidual_ -= localResidual_.residual();
            derivStorage_[scvIdx] -= localResidual_.storageTerm()[scvIdx];
        }
        else {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            derivResidual_ -= residual_;
            derivStorage_[scvIdx] -= residualStorage_[scvIdx];
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        derivResidual_ /= delta;
        derivStorage_[scvIdx] /= delta;

        // restore the original state of the element's volume
        // variables
        elemVars.restoreScvVars(scvIdx);

#ifndef NDEBUG
        for (unsigned i = 0; i < derivResidual_.size(); ++i)
            Valgrind::CheckDefined(derivResidual_[i]);
#endif
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at vertex 'scvIdx' .
     */
    void updateLocalJacobian_(const ElementVariables &elemVars,
                              int scvIdx,
                              int pvIdx)
    {
        // store the derivative of the storage term
        for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
            jacobianStorage_[scvIdx][eqIdx][pvIdx] = derivStorage_[scvIdx][eqIdx];
        }

        int numScv = elemVars.numScv();
        for (int eqScvIdx = 0; eqScvIdx < numScv; eqScvIdx++)
        {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[eqScvIdx][scvIdx][eqIdx][pvIdx] is the rate of
                // change of the residual of equation 'eqIdx' at
                // vertex 'eqScvIdx' depending on the primary variable
                // 'pvIdx' at vertex 'scvIdx'.
                jacobian_[eqScvIdx][scvIdx][eqIdx][pvIdx] = derivResidual_[eqScvIdx][eqIdx];
                Valgrind::CheckDefined(jacobian_[eqScvIdx][scvIdx][eqIdx][pvIdx]);
            }
        }
    }

    Problem *problemPtr_;
    Model *modelPtr_;

    ElementVariables *internalElemVars_;

    LocalBlockMatrix jacobian_;
    LocalStorageMatrix jacobianStorage_;

    LocalBlockVector residual_;
    LocalBlockVector residualStorage_;

    LocalBlockVector derivResidual_;
    LocalBlockVector derivStorage_;

    LocalResidual localResidual_;
};
}

#endif
