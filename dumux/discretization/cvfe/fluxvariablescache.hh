// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Flux variables cache class for control-volume finite element schemes
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_FLUXVARIABLES_CACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_FLUXVARIABLES_CACHE_HH

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variables cache class for control-volume finite element schemes.
 *        For control-volume finite element schemes, this class does not contain any physics-/process-dependent
 *        data. It solely stores disretization-/grid-related data.
 */
template< class Scalar, class GridGeometry >
class CVFEFluxVariablesCache
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;

public:
    //! whether the cache needs an update when the solution changes
    static bool constexpr isSolDependent = false;

    //! update the cache for an scvf
    template< class Problem, class ElementVolumeVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        update(problem, element, fvGeometry, elemVolVars, scvf.ipGlobal());
    }

    //! update the cache for a given global position
    template< class Problem, class ElementVolumeVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const GlobalPosition& globalPos)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        ipGlobal_ = globalPos;
        const auto ipLocal = geometry.local(ipGlobal_);
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal, shapeValues_); // shape values for rho

        // compute the gradN at for every scv/dof
        gradN_.resize(fvGeometry.numScv());
        for (const auto& scv: scvs(fvGeometry))
            jacInvT_.mv(shapeJacobian_[scv.localDofIndex()][0], gradN_[scv.indexInElement()]);
    }

    //! returns the global position for which this cache has been updated
    const GlobalPosition& ipGlobal() const { return ipGlobal_; }
    //! returns the shape function gradients in local coordinates at the integration point
    const std::vector<ShapeJacobian>& shapeJacobian() const { return shapeJacobian_; }
    //! returns the shape function values at the integration point
    const std::vector<ShapeValue>& shapeValues() const { return shapeValues_; }
    //! returns inverse transposed jacobian at the integration point
    const JacobianInverseTransposed& jacInvT() const { return jacInvT_; }
    //! returns the shape function gradients in global coordinates at the integration point
    const GlobalPosition& gradN(unsigned int scvIdxInElement) const { return gradN_[scvIdxInElement]; }

private:
    GlobalPosition ipGlobal_;
    std::vector<GlobalPosition> gradN_;
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;
};

} // end namespace Dumux

#endif
