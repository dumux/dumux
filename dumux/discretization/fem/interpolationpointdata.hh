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
    using ShapeGradients = std::vector<GlobalPosition>;

public:
    // The default constructor
    FEInterpolationPointData() = delete;

    // The constructor
    template<class Geometry>
    FEInterpolationPointData(const Geometry& geometry,
                             const LocalPosition& local,
                             const LocalBasis& localBasis)
    : local_(local),
      global_(geometry.global(local))
    {
        auto numLocalDofs = localBasis.size();

        // set the shape values
        shapeValues_.resize(numLocalDofs);
        localBasis.evaluateFunction(local, shapeValues_);

        // the local shape function gradients
        std::vector<JacobianType> shapeGrads(numLocalDofs);
        localBasis.evaluateJacobian(local, shapeGrads);

        // the global shape function gradients
        const auto jacInvT = geometry.jacobianInverseTransposed(local);
        shapeGradients_.resize(numLocalDofs, GlobalPosition(0.0));
        for (unsigned int i = 0; i < numLocalDofs; ++i)
            jacInvT.umv(shapeGrads[i][0], shapeGradients_[i]);
    }

    //! The shape values at the quadrature point
    const ShapeValues& shapeValues() const
    { return shapeValues_; }

    //! The shape value of a local dof at the quadrature point
    const RangeType& shapeValue(int i) const
    { return shapeValues_[i]; }

    //! The shape value gradients at the quadrature point
    const ShapeGradients& shapeGradients() const
    { return shapeGradients_; }

    //! The shape value gradient of a local dof at the quadrature point
    const GlobalPosition& gradN(int i) const
    { return shapeGradients_[i]; }

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
    ShapeGradients shapeGradients_;
};

/*!
 * \ingroup FEMDiscretization
 * \brief Interpolation point data related to a face of an element
 */
template<class GlobalPosition, class LocalBasis, class BoundaryFlag>
class FEFaceInterpolationPointData : public FEInterpolationPointData<GlobalPosition, LocalBasis>
{
    using ParentType = FEInterpolationPointData<GlobalPosition, LocalBasis>;
    using LocalPosition = typename LocalBasis::Traits::DomainType;
public:
    // The default constructor
    FEFaceInterpolationPointData() = delete;

    // The constructor
    template<class Geometry>
    FEFaceInterpolationPointData(const Geometry& geometry,
                                 const LocalPosition& local,
                                 const LocalBasis& localBasis,
                                 const GlobalPosition& n,
                                 const BoundaryFlag& bFlag)
    : ParentType(geometry, local, localBasis), normal_(n), boundaryFlag_(bFlag)
    {}

    //! The unit outer normal vector at the quadrature point
    const GlobalPosition& unitOuterNormal() const
    { return normal_; }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

private:
    const GlobalPosition& normal_;
    BoundaryFlag boundaryFlag_;
};


} // end namespace Dumux

#endif
