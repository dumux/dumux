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
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using JacobianAssembler = typename GET_PROP_TYPE(TypeTag, JacobianAssembler);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Model = typename GET_PROP_TYPE(TypeTag, Model);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    // copying a local jacobian is not a good idea
    ImplicitLocalJacobian(const ImplicitLocalJacobian &) = delete;

    ImplicitLocalJacobian()
    {
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
    void assemble(const Element& element, JacobianMatrix& matrix)
    {
        DUNE_THROW(Dune::NotImplemented, "Assemble routine not provided by the actual implementation of the local jacobian!");
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
     * \brief Returns the residual of the equations at dof i.
     *
     * \param i The local dof index on which
     *          the equations are defined
     */
    const PrimaryVariables &residual(const int i) const
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
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the model.
     */
    Model &model_()
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


    Scalar numericEpsilon(const Scalar priVar) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
    }

    /*!
     * \brief Updates the current global Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col'.
     */
    void updateGlobalJacobian_(JacobianMatrix& matrix,
                               const int globalI,
                               const int globalJ,
                               const int pvIdx,
                               const PrimaryVariables &partialDeriv)
    {
//	std::cout << std::endl;
//	printmatrix(std::cout, matrix, "nachAsssembleMatrix", "", 7, 0);
        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
        {
            // A[i][col][eqIdx][pvIdx] is the rate of change of
            // the residual of equation 'eqIdx' at dof 'i'
            // depending on the primary variable 'pvIdx' at dof
            // 'col'.
            matrix[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
            Valgrind::CheckDefined(matrix[globalI][globalJ][eqIdx][pvIdx]);
        }
    }

    // The problem we would like to solve
    Problem *problemPtr_;

    LocalResidual localResidual_;

    ElementSolutionVector residual_;
};

}

#endif
