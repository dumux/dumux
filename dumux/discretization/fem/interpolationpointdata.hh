// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEMDiscretization
 * \brief Shape functions and gradients at an interpolation point
 */
#ifndef DUMUX_DISCRETIZATION_FE_INTEGRATION_POINT_DATA_HH
#define DUMUX_DISCRETIZATION_FE_INTEGRATION_POINT_DATA_HH

#include <vector>

namespace Dumux {

template<class GlobalPosition, class LocalBasis>
class FEInterpolationPointData
{
    using LocalPosition = typename LocalBasis::Traits::DomainType;
    using RangeType = typename LocalBasis::Traits::RangeType;
    using JacobianType = typename LocalBasis::Traits::JacobianType;

    using ShapeValues = std::vector<RangeType>;
    using ShapeJacobians = std::vector<JacobianType>;
    using Gradients = std::vector<GlobalPosition>;

public:
    // The default constructor
    FEInterpolationPointData() = delete;

    // The constructor
    template<class Geometry>
    FEInterpolationPointData(const Geometry& geometry,
                             const LocalPosition& local,
                             const GlobalPosition& global,
                             const LocalBasis& localBasis)
    : local_(local),
      global_(global)
    {
        auto numLocalDofs = localBasis.size();

        // set the shape values
        shapeValues_.resize(numLocalDofs);
        localBasis.evaluateFunction(local, shapeValues_);

        // the local shape function gradients
        shapeJacobians_.resize(numLocalDofs);
        localBasis.evaluateJacobian(local, shapeJacobians_);

        // the global shape function gradients
        const auto jacInvT = geometry.jacobianInverseTransposed(local);
        gradients_.resize(numLocalDofs, GlobalPosition(0.0));
        for (unsigned int i = 0; i < numLocalDofs; ++i)
            jacInvT.umv(shapeJacobians_[i][0], gradients_[i]);
    }

    // The constructor
    template<class Geometry>
    FEInterpolationPointData(const Geometry& geometry,
                             const LocalPosition& local,
                             const LocalBasis& localBasis)
    : FEInterpolationPointData(geometry, local, geometry.global(local), localBasis)
    { }

    //! The shape values at the quadrature point
    const ShapeValues& shapeValues() const
    { return shapeValues_; }

    //! The shape value of a local dof at the quadrature point
    const RangeType& shapeValue(int i) const
    { return shapeValues_[i]; }

    //! The shape value gradients at the quadrature point
    const ShapeJacobians& shapeJacobians() const
    { return shapeJacobians_; }

    //! The shape value gradient of a local dof at the quadrature point
    const GlobalPosition& gradN(int i) const
    { return gradients_[i]; }

    //! The local position of the quadrature point
    const LocalPosition& local() const
    { return local_; }

    //! The global position of the quadrature point
    const GlobalPosition& global() const
    { return global_; }

private:
    LocalPosition local_;
    GlobalPosition global_;

    ShapeValues shapeValues_;
    ShapeJacobians shapeJacobians_;
    Gradients gradients_;
};

/*!
 * \ingroup FEMDiscretization
 * \brief An interpolation point related to an intersection
 */
template<class BaseClass, class BoundaryFlag, class LocalIndex>
class FEFaceInterpolationPointData : public BaseClass
{
public:
    using GlobalPosition = std::remove_cvref_t<decltype(std::declval<BaseClass>().global())>;
    using LocalPosition = std::remove_cvref_t<decltype(std::declval<BaseClass>().local())>;

    template<class... Args>
    FEFaceInterpolationPointData(GlobalPosition&& n, BoundaryFlag&& bFlag, LocalIndex index, Args&&... args)
    : BaseClass(std::forward<Args>(args)...), normal_(std::move(n)), boundaryFlag_(std::move(bFlag)), index_(index) {}

    template<class... Args>
    FEFaceInterpolationPointData(const GlobalPosition& n, const BoundaryFlag& bFlag, LocalIndex index, Args&&... args)
    : BaseClass(std::forward<Args>(args)...), normal_(n), boundaryFlag_(bFlag), index_(index) {}

    //! The unit outer normal vector at the quadrature point
    const GlobalPosition& unitOuterNormal() const
    { return normal_; }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

    //! The local index of an intersection (default is intersection.indexInInside())
    LocalIndex intersectionIndex() const
    { return index_; }

private:
    const GlobalPosition& normal_;
    BoundaryFlag boundaryFlag_;
    LocalIndex index_;
};

} // end namespace Dumux

#endif
