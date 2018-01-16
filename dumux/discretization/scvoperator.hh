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
 *
 * \ingroup Discretization
 * \brief Class for the evaluation of primary variables and gradients on sub-control volumes.
 */
#ifndef DUMUX_SUBCONTROL_VOLUME_OPERATOR_HH
#define DUMUX_SUBCONTROL_VOLUME_OPERATOR_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
//! Forward declaration of the method-specific implementation
template<class TypeTag, bool isBox>
class SubControlVolumeOperatorImplementation;

/**
 * \brief The base class for solution dependent spatial parameters.
 */
template<class TypeTag>
using SubControlVolumeOperator DUNE_DEPRECATED_MSG("Use evalSolution() instead") =
        SubControlVolumeOperatorImplementation<TypeTag, (GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box)>;


//! Specialization for the box method
template<class TypeTag>
class SubControlVolumeOperatorImplementation<TypeTag, true>
{
    using ThisType = SubControlVolumeOperatorImplementation<TypeTag, true>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using CoordScalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using PriVarGradients = Dune::FieldVector<GlobalPosition, numEq>;

public:

    /*!
     * \brief Interpolates a given solution to the scv center.
     *
     * \return the interpolated Primary Variables
     * \param element The element
     * \param scv The sub-control volume
     * \param elementSolution The primary variables at the dofs of the element.
     */
    static PrimaryVariables evaluateSolution(const Element& element,
                                             const SubControlVolume& scv,
                                             const ElementSolutionVector& elementSolution)
    {
        auto g = element.geometry();

        // The local basis of the finite element
        FeCache feCache;
        const auto& localBasis = feCache.get(g.type()).localBasis();

        // evaluate the shape functions at the scv center
        auto local = g.local(scv.center());
        std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;
        localBasis.evaluateFunction(local, shapeValues);

        // interpolate the solution
        PrimaryVariables result(0.0);
        for (int i = 0; i < element.subEntities(dim); ++i)
        {
            auto value = elementSolution[i];
            value *= shapeValues[i][0];
            result += value;
        }

        return result;
    }

    /*!
     * \brief Evaluates the gradient of a given solution at the scv center.
     *
     * \return the interpolated Primary Variables
     * \param element The element
     * \param scv The sub-control volume
     * \param elementSolution The primary variables at the dofs of the element.
     */
    static PriVarGradients evaluateGradients(const Element& element,
                                             const SubControlVolume& scv,
                                             const ElementSolutionVector& elementSolution)
    {
        auto g = element.geometry();

        // The local basis of the finite element
        FeCache feCache;
        const auto& localBasis = feCache.get(g.type()).localBasis();

        // the scv center in local coordinates
        auto local = g.local(scv.center());

        // the inverse transposed of the jacobian matrix
        const auto jacInvT = g.jacobianInverseTransposed(local);

        // evaluate the shape function gradients at the local scv center
        std::vector< ShapeJacobian > shapeJacobian;
        localBasis.evaluateJacobian(local, shapeJacobian);

        // interpolate the gradients
        PriVarGradients result( PrimaryVariables(0.0) );
        for (int i = 0; i < element.subEntities(dim); ++i)
        {
            // the global shape function gradient
            GlobalPosition gradN;
            jacInvT.mv(shapeJacobian[i][0], gradN);

            // add gradient to global privar gradients
            for (unsigned int pvIdx = 0; pvIdx < numEq; ++pvIdx)
            {
                GlobalPosition tmp(gradN);
                tmp *= elementSolution[i][pvIdx];
                result[pvIdx] += tmp;
            }
        }

        return result;
    }
};

// Specialization for cell-centered schemes
template<class TypeTag>
class SubControlVolumeOperatorImplementation<TypeTag, false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    using Element = typename GridView::template Codim<0>:: Entity;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using PriVarGradients = Dune::FieldVector<GlobalPosition, numEq>;

public:

    /*!
     * \brief Evaluates a given solution at the scv center
     *
     * \return the interpolated Primary Variables
     * \param element The element
     * \param scv The sub-control volume
     * \param elementSolution The primary variables at the dofs of the element.
     */
    static PrimaryVariables evaluateSolution(const Element& element,
                                             const SubControlVolume& scv,
                                             const ElementSolutionVector& elementSolution)
    { return elementSolution[0]; }

    /*!
     * \brief Evaluates the gradient of a given solution at the scv center.
     *
     * \return the interpolated Primary Variables
     * \param element The element
     * \param scv The sub-control volume
     * \param elementSolution The primary variables at the dofs of the element.
     *
     * TODO how could this be achieved for cc methods?
     */
    static PrimaryVariables evaluateGradients(const Element& element,
                                              const SubControlVolume& scv,
                                              const ElementSolutionVector& elementSolution)
    { DUNE_THROW(Dune::NotImplemented, "Gradient evaluation at scv centers for cell-centered methods"); }
};

} // namespace Dumux

#endif
