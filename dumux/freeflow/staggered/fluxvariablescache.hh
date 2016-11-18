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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREEFLOW_IMPLICIT_FLUXVARIABLESCACHE_HH
#define DUMUX_FREEFLOW_IMPLICIT_FLUXVARIABLESCACHE_HH

#include <dumux/implicit/properties.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethods Method>
class FreeFlowFluxVariablesCacheImplementation
{};

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store flux stencils and data required for flux calculation
 */
template<class TypeTag>
using FreeFlowFluxVariablesCache = FreeFlowFluxVariablesCacheImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

// specialization for the Box Method
template<class TypeTag>
class FreeFlowFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;
    using TransmissibilityVector = std::vector<IndexType>;

    using CoordScalar = typename GridView::ctype;
    static const int dim = GridView::dimension;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;

public:

    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const SubControlVolumeFace &scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvf.center());
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        //localBasis.evaluateFunction(ipLocal, shapeValue_); // do we need the shapeValues for the flux?

        // The stencil info is obsolete for the box method.
        // It is here for compatibility with cc methods
        stencil_ = Stencil(0);
    }

    const std::vector<ShapeJacobian>& shapeJacobian() const
    { return shapeJacobian_; }

   /* const std::vector<ShapeValue>& shapeValue() const
    { return shapeValue_; }*/

    const JacobianInverseTransposed& jacInvT() const
    { return jacInvT_; }

    const Stencil& stencil() const
    {
        return stencil_;
    }

private:
    std::vector<ShapeJacobian> shapeJacobian_;
    //std::vector<ShapeValue> shapeValue_;
    JacobianInverseTransposed jacInvT_;

    Stencil stencil_;
};

// specialization for the cell centered tpfa method
template<class TypeTag>
class FreeFlowFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeStencil(problem, element, fvGeometry, scvf);
        tij_ = AdvectionType::calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);
    }

    const Stencil& stencil() const
    { return stencil_; }

    const Scalar& tij() const
    { return tij_; }

private:
    Stencil stencil_;
    Scalar tij_;
};

// specialization for the cell centered tpfa method
template<class TypeTag>
class FreeFlowFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::Staggered>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeCellCenterStencil(problem, element, fvGeometry, scvf);
        tij_ = AdvectionType::calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);
    }

    const Stencil& stencil() const
    { return stencil_; }

    const Scalar& tij() const
    { return tij_; }

private:
    Stencil stencil_;
    Scalar tij_;
};

} // end namespace

#endif
