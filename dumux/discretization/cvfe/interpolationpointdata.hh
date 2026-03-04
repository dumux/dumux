// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Classes representing interpolation point data for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_IP_DATA_HH
#define DUMUX_CVFE_IP_DATA_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to an element that includes global and local positions
 */
template<class LocalPos, class GlobalPos>
class InterpolationPointData
{
public:
    using LocalPosition = LocalPos;
    using GlobalPosition = GlobalPos;

    InterpolationPointData(LocalPosition&& localPos, GlobalPosition&& pos)
    : local_(std::move(localPos)), global_(std::move(pos)) {}
    InterpolationPointData(const LocalPosition& localPos, const GlobalPosition& pos)
    :  local_(localPos), global_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& global() const
    { return global_; }

    //! The local position of the quadrature point
    const LocalPosition& local() const
    { return local_; }

private:
    LocalPosition local_;
    GlobalPosition global_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to a localDof of an element, giving its global and local positions
 */
template<class LocalPosition, class GlobalPosition, class LocalIndex>
class LocalDofInterpolationPointData : public InterpolationPointData<LocalPosition, GlobalPosition>
{
    using ParentType = InterpolationPointData<LocalPosition, GlobalPosition>;
public:
    LocalDofInterpolationPointData(LocalPosition&& localPos, GlobalPosition&& pos, LocalIndex index)
    : ParentType(localPos, pos), localDofIndex_(index) {}
    LocalDofInterpolationPointData(const LocalPosition& localPos, const GlobalPosition& pos, LocalIndex index)
    : ParentType(localPos, pos), localDofIndex_(index) {}

    //! The local index of the corresponding dof
    LocalIndex localDofIndex() const
    { return localDofIndex_; }

private:
    LocalIndex localDofIndex_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to a global position of an element, giving its local positions by a mapping
 */
template<class LocalMapping, class GlobalPos>
class InterpolationPointDataLocalMapping
{
public:
    using LocalPosition = std::invoke_result_t<LocalMapping, const GlobalPos&>;
    using GlobalPosition = GlobalPos;

    InterpolationPointDataLocalMapping(LocalMapping&& mapping, GlobalPosition&& pos) : localMapping_(std::move(mapping)), global_(std::move(pos)) {}
    InterpolationPointDataLocalMapping(LocalMapping&& mapping, const GlobalPosition& pos) : localMapping_(std::move(mapping)), global_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& global() const
    { return global_; }

    //! The local position of the quadrature point
    const LocalPosition local() const
    { return localMapping_(global_); }

private:
    LocalMapping localMapping_;
    GlobalPosition global_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to a face of an element
 */
template<class BaseClass, class LocalIndex>
class FaceInterpolationPointData : public BaseClass
{
public:
    using GlobalPosition = std::remove_cvref_t<decltype(std::declval<BaseClass>().global())>;
    using LocalPosition = std::remove_cvref_t<decltype(std::declval<BaseClass>().local())>;

    template<class... Args>
    FaceInterpolationPointData(GlobalPosition&& n, LocalIndex index, Args&&... args)
    : BaseClass(std::forward<Args>(args)...), normal_(std::move(n)), scvfIndex_(index) {}

    template<class... Args>
    FaceInterpolationPointData(const GlobalPosition& n, LocalIndex index, Args&&... args)
    : BaseClass(std::forward<Args>(args)...), normal_(n), scvfIndex_(index) {}

    //! The unit outer normal vector at the quadrature point
    const GlobalPosition& unitOuterNormal() const
    { return normal_; }

    //! The local index of an scvf
    LocalIndex scvfIndex() const
    { return scvfIndex_; }

private:
    GlobalPosition normal_;
    LocalIndex scvfIndex_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Wraps interpolation point data and adds a quadrature point index for use in quadrature loops
 */
template<class IpData>
class IndexedQuadratureInterpolationPointData : public IpData
{
public:
    template<class... Args>
    IndexedQuadratureInterpolationPointData(std::size_t qpIdx, Args&&... args)
    : IpData(std::forward<Args>(args)...), qpIndex_(qpIdx) {}

    //! The quadrature point index
    std::size_t qpIndex() const
    { return qpIndex_; }

private:
    std::size_t qpIndex_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Interpolation point data related to a local basis
 */
template< class GridGeometry >
class LocalBasisInterpolationPointData
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using LocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using RangeType = typename LocalBasis::Traits::RangeType;
    using JacobianType = typename LocalBasis::Traits::JacobianType;

    using Gradients = std::vector<GlobalPosition>;

    using FVElementGeometry = typename GridGeometry::LocalView;

public:
    //! whether the cache needs an update when the solution changes
    static constexpr bool isSolDependent = false;

    //! update the cache for a given global position
    template< class Problem, class ElementVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVariables& elemVars,
                const GlobalPosition& globalPos)
    {
        update_(fvGeometry, fvGeometry.elementGeometry().local(globalPos));
    }

    //! update the cache for interpolation point data
    template< class Problem, class ElementVariables, Concept::IpData IpData >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVariables& elemVars,
                const IpData& ipData)
    {
        update_(fvGeometry, ipData.local());
    }

    //! returns the shape function gradients in local coordinates at the interpolation point
    const std::vector<JacobianType>& shapeJacobian() const { return shapeJacobian_; }
    //! returns the shape function values at the interpolation point
    const std::vector<RangeType>& shapeValues() const { return shapeValues_; }
    //! returns the shape function gradients in global coordinates at the interpolation point
    const GlobalPosition& gradN(unsigned int localIdx) const { return gradN_[localIdx]; }

private:
    //! update the cache for a given local and global position
    template<class LocalPosition>
    void update_(const FVElementGeometry& fvGeometry,
                 const LocalPosition& localPos)
    {
        const auto& geometry = fvGeometry.elementGeometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the interpolation point
        const auto jacInvT = geometry.jacobianInverseTransposed(localPos);
        localBasis.evaluateJacobian(localPos, shapeJacobian_);
        localBasis.evaluateFunction(localPos, shapeValues_); // shape values for rho

        // compute the gradN for every local dof
        gradN_.resize(Detail::LocalDofs::numLocalDofs(fvGeometry));
        for (const auto& localDof: localDofs(fvGeometry))
            jacInvT.mv(shapeJacobian_[localDof.index()][0], gradN_[localDof.index()]);
    }

    Gradients gradN_;
    std::vector<JacobianType> shapeJacobian_;
    std::vector<RangeType> shapeValues_;
};

} // end namespace Dumux::CVFE

#endif
