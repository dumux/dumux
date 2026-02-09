// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Flux variables cache class for control-volume finite element schemes
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_FLUXVARIABLES_CACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_FLUXVARIABLES_CACHE_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

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
    using LocalPosition = typename Element::Geometry::LocalCoordinate;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;

    using IpData = Dumux::CVFE::InterpolationPointData<LocalPosition, GlobalPosition>;

public:
    using ScvfQuadratureRule = CVFE::Detail::ScvfQuadratureRuleType<FVElementGeometry>;

    //! whether the cache needs an update when the solution changes
    static bool constexpr isSolDependent = false;

    //! update the cache for an scvf (only enabled for MidpointQuadrature)
    template< class Problem, class ElementVolumeVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        if constexpr (std::is_same_v<ScvfQuadratureRule, QuadratureRules::MidpointQuadrature>)
            this->update(problem, element, fvGeometry, elemVolVars, scvf.ipGlobal());
        else
            DUNE_THROW(Dune::NotImplemented, "CVFEFluxVariablesCache: update for scvf only implemented for MidpointQuadrature.");
    }

    //! update the cache for a given global position
    template< class Problem, class ElementVolumeVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const GlobalPosition& globalPos)
    {
        const auto& geometry = elementGeometry_(fvGeometry);
        ipGlobal_ = globalPos;
        ipLocal_ = geometry.local(globalPos);
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal_);
        evaluateBasis_(fvGeometry);
    }

    //! update the cache for interpolation point data
    template< class Problem, class ElementVolumeVariables, Concept::IpData IpData >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const IpData& ipData)
    {
        const auto& geometry = elementGeometry_(fvGeometry);
        ipGlobal_ = ipData.global();
        ipLocal_ = ipData.local();
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal_);
        evaluateBasis_(fvGeometry);
    }

    //! returns the global position for which this cache has been updated
    const GlobalPosition& ipGlobal() const { return ipGlobal_; }
    //! returns the local position for which this cache has been updated
    const LocalPosition& ipLocal() const { return ipLocal_; }
    //! returns the ipData for which this cache has been updated
    IpData ipData() const { return IpData(ipLocal(), ipGlobal()); }
    //! returns the shape function gradients in local coordinates at the interpolation point
    const std::vector<ShapeJacobian>& shapeJacobian() const { return shapeJacobian_; }
    //! returns the shape function values at the interpolation point
    const std::vector<ShapeValue>& shapeValues() const { return shapeValues_; }
    //! returns inverse transposed jacobian at the interpolation point
    const JacobianInverseTransposed& jacInvT() const { return jacInvT_; }
    //! returns the shape function gradients in global coordinates at the interpolation point
    const GlobalPosition& gradN(unsigned int scvIdxInElement) const { return gradN_[scvIdxInElement]; }

private:
    static decltype(auto) elementGeometry_(const FVElementGeometry& fvGeometry)
    {
        if constexpr (requires { fvGeometry.elementGeometry(); })
            return fvGeometry.elementGeometry();
        else
            return fvGeometry.element().geometry();
    }

    //! update the cache for a given local and global position
    void evaluateBasis_(const FVElementGeometry& fvGeometry)
    {
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the interpolation point
        localBasis.evaluateJacobian(ipLocal_, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal_, shapeValues_); // shape values for rho

        // compute the gradN at for every scv/dof
        gradN_.resize(Detail::LocalDofs::numLocalDofs(fvGeometry));
        for (const auto& localDof: localDofs(fvGeometry))
            jacInvT_.mv(shapeJacobian_[localDof.index()][0], gradN_[localDof.index()]);
    }

    std::vector<GlobalPosition> gradN_;
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;
    GlobalPosition ipGlobal_;
    LocalPosition ipLocal_;
};

} // end namespace Dumux

#endif
