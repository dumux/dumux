// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Caculates the Jacobian of the local residual for fully-implicit models
 */
#ifndef DUMUX_IMPLICIT_LOCAL_JACOBIAN_HH
#define DUMUX_IMPLICIT_LOCAL_JACOBIAN_HH

#include <dune/istl/io.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ImplicitLocalJacobian
 * \brief Calculates the Jacobian of the local residual for fully-implicit models
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
class ImplicitLocalJacobian
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension,

        Green = JacobianAssembler::Green
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    // copying a local jacobian is not a good idea
    ImplicitLocalJacobian(const ImplicitLocalJacobian &);

public:
    ImplicitLocalJacobian()
    {
        numericDifferenceMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);
        Valgrind::SetUndefined(problemPtr_);
    }


    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;
        localResidual_.init(problem);
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
        reset_();

        bcTypes_.update(problem_(), element_(), fvElemGeom_());

        // calculate the local residual
        localResidual().eval(element_(), bcTypes_);
        residual_ = localResidual().residual();

        model_().updatePVWeights(fvElemGeom_());

        // get stencil informations
        const auto& elementStencil = model_().stencils().elementStencil(element);

        // set size of local jacobian matrix
        const std::size_t numCols = elementStencil.size();
        std::size_t numRows;

        if (isBox)
            numRows = numCols;
        else
            numRows = 1;

        A_.setSize(numRows, numCols);

        // calculate derivatives for the dofs inside the element
        ElementSolutionVector partialDeriv(numRows);
        for (auto&& scv : fvElemGeom_().scvs())
        {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                asImp_().evalPartialDerivative_(partialDeriv, scv, pvIdx);

                // update the local stiffness matrix with the current partial derivatives
                updateLocalJacobian_(scv.indexInElement(), pvIdx, partialDeriv);
            }
        }

        // TODO: calculate derivatives in the case of an extended source stencil
        // const auto& extendedSourceStencil = model_().stencils().extendedSourceStencil(element);
        // for (auto&& globalJ : extendedSourceStencil)
        // {
        //     for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
        //         {
        //             evalPartialDerivativeSource_(partialDeriv, globalJ, pvIdx, neighborToFluxVars[globalJ]);

        //             // update the local stiffness matrix with the partial derivatives
        //             updateLocalJacobian_(j++, pvIdx, partialDeriv);
        //         }
        // }

        // for cellcentered methods, calculate the derivatives w.r.t cells in stencil
        if (!isBox)
        {
            const auto& neighborStencil = model_().stencils().neighborStencil(element);

            // map each neighbor dof to a set of fluxVars that need to be recalculated
            // i.o.t calculate derivative w.r.t this neighbor
            std::map< unsigned int, std::set<unsigned int> > neighborToFluxVars;

            // loop over scvFaces/fluxVars of the element
            for (auto&& scvFace : fvElemGeom_().scvf())
            {
                int fluxVarIdx = scvFace.index();
                const auto& fluxVars = model_().fluxVars(fluxVarIdx);
                for (auto&& globalJ : neighborStencil)
                    if (fluxVars.stencil().count(globalJ))
                        neighborToFluxVars[globalJ].insert(fluxVarIdx);
            }

            // loop over the neighbors and calculation of the change in flux
            // with a change in the primary variables at the neighboring dof
            int j = 1;
            for (auto&& globalJ : neighborStencil)
            {
                for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
                {
                    evalPartialDerivativeFlux_(partialDeriv, globalJ, pvIdx, neighborToFluxVars[globalJ]);

                    // update the local stiffness matrix with the partial derivatives
                    updateLocalJacobian_(j++, pvIdx, partialDeriv);
                }
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
     * \brief Returns the Jacobian of the equations at subcontrolvolume i
     * to the primary variables at subcontrolvolume j.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     * \param j The local subcontrolvolume index which holds
     *          primary variables
     */
    const MatrixBlock &mat(const int i, const int j) const
    { return A_[i][j]; }

    /*!
     * \brief Returns the residual of the equations at subcontrolvolume i.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     */
    const PrimaryVariables &residual(const int i) const
    { return residual_[i]; }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param scvIdx     The local index of the element's subcontrolvolume for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon(const SubControlVolume &scv,
                          const int pvIdx) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        Scalar priVar = model_().curVolVars(scv).priVar(pvIdx);
        return baseEps*(std::abs(priVar) + 1.0);
    }

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
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem_()
    {
        Valgrind::CheckDefined(problemPtr_);
        return *problemPtr_;
    }

    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the element.
     */
    const Element &element_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the jacobian assembler.
     */
    const JacobianAssembler &jacAsm_() const
    { return model_().jacobianAssembler(); }

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Reset the local jacobian matrix to 0
     */
    void reset_()
    { A_ = 0.0; }

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
     * \param partialDeriv The vector storing the partial derivatives of all
     *              equations
     * \param col The block column index of the degree of freedom
     *            for which the partial derivative is calculated.
     *            Box: a sub-control volume index.
     *            Cell centered: a neighbor index.
     * \param pvIdx The index of the primary variable
     *              for which the partial derivative is calculated
     */
    void evalPartialDerivative_(ElementSolutionVector &partialDeriv,
                                const SubControlVolume &scv,
                                const int pvIdx)
    {
        int dofIdxGlobal = scv.dofIndex();

        auto priVars = model_().curSol()[dofIdxGlobal];
        auto origVolVars = model_().curVolVars(scv);

        Scalar eps = asImp_().numericEpsilon(scv, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // update the volume variables
            model_().curVolVars(scv).update(priVars, problem_(), element_(), scv);

            // calculate the residual with the deflected primary variables
            localResidual().eval(element_(), bcTypes_);

            // store the residual and the storage term
            partialDeriv = localResidual().residual();
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = residual_;
        }

        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // update the volume variables
            model_().curVolVars(scv).update(priVars, problem_(), element_(), scv);

            // calculate the residual with the deflected primary variables
            localResidual().eval(element_(), bcTypes_);

            // subtract the residual from the derivative storage
            partialDeriv -= localResidual().residual();
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= residual_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

        // restore the original state of the scv's volume variables
        model_().curVolVars_(scv) = origVolVars;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }

    void evalPartialDerivativeFlux_(ElementSolutionVector &partialDeriv,
                                    const unsigned int globalJ,
                                    const int pvIdx,
                                    const std::set<unsigned int> &fluxVarsJ)
    {
        if (isBox)
            DUNE_THROW(Dune::InvalidStateException, "Calling evalPartialDerivativeFlux_(...) for box method.");

        auto&& scvJ = model_().fvGeometries().subControlVolume(globalJ);
        auto priVarsJ = model_().curSol()[globalJ];
        auto origVolVarsJ = model_().curVolVars(scvJ);

        // calculate the flux in the undeflected state
        PrimaryVariables origFlux = 0.0;
        for (auto&& fluxVarIdx : fluxVarsJ)
            origFlux += localResidual().evalFlux_(fluxVarIdx);

        Scalar eps = asImp_().numericEpsilon(scvJ, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_ >= 0)
        {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVarsJ[pvIdx] += eps;
            delta += eps;

            // update the volume variables
            model_().curVolVars(scvJ).update(priVarsJ, problem_(), element_(), scvJ);

            // calculate the flux with the deflected primary variables
            // TODO: for solution dependent spatial params fluxVar update needed!
            PrimaryVariables deflectFlux = 0.0;
            for (auto&& fluxVarIdx : fluxVarsJ)
                deflectFlux += localResidual().evalFlux_(fluxVarIdx);

            // store the calculated flux
            partialDeriv = deflectFlux;
        }
        else
        {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) flux f(x)
            partialDeriv = PrimaryVariables;
        }

        if (numericDifferenceMethod_ <= 0)
        {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVarsJ[pvIdx] -= delta + eps;
            delta += eps;

            // update the volume variables
            model_().curVolVars(scvJ).update(priVarsJ, problem_(), element_(), scvJ);

            // calculate the flux with the deflected primary variables
            // TODO: for solution dependent spatial params fluxVar update needed!
            PrimaryVariables deflectFlux = 0.0;
            for (auto&& fluxVarIdx : fluxVarsJ)
                PrimaryVariables += localResidual().evalFlux_(fluxVarIdx);

            // subtract the residual from the derivative storage
            partialDeriv -= deflectFlux;
        }
        else
        {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) flux f(x)
            partialDeriv -= origFlux;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;

        // restore the original state of the scv's volume variables
        model_().curVolVars_(scvJ) = origVolVarsJ;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col' .
     */
    void updateLocalJacobian_(const int col,
                              const int pvIdx,
                              const ElementSolutionVector &partialDeriv)
    {
        for (auto&& scv : fvElemGeom_().scvs())
        {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                int i = scv.indexInElement();
                // A[i][col][eqIdx][pvIdx] is the rate of change of
                // the residual of equation 'eqIdx' at dof 'i'
                // depending on the primary variable 'pvIdx' at dof
                // 'col'.
                this->A_[i][col][eqIdx][pvIdx] = partialDeriv[i][eqIdx];
                Valgrind::CheckDefined(this->A_[i][col][eqIdx][pvIdx]);
            }
        }
    }

    const FVElementGeometry& fvElemGeom_() const
    {
        return model_().fvGeometries(problem_().elementMapper().index(element_()));
    }

    const Element *elemPtr_;
    ElementBoundaryTypes bcTypes_;

    // The problem we would like to solve
    Problem *problemPtr_;

    LocalResidual localResidual_;

    LocalBlockMatrix A_;

    ElementSolutionVector residual_;

    int numericDifferenceMethod_;
};
}

#endif
