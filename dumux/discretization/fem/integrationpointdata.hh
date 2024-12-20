// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEMDiscretization
 * \brief Shape functions and gradients at an integration point
 */
#ifndef DUMUX_DISCRETIZATION_FE_INTEGRATION_POINT_DATA_HH
#define DUMUX_DISCRETIZATION_FE_INTEGRATION_POINT_DATA_HH

#include <vector>

namespace Dumux {

template<class GlobalPosition, class LocalBasis>
class FEIntegrationPointData
{
    using LocalPosition = typename LocalBasis::Traits::DomainType;
    using RangeType = typename LocalBasis::Traits::RangeType;
    using JacobianType = typename LocalBasis::Traits::JacobianType;

    using ShapeValues = std::vector<RangeType>;
    using ShapeGradients = std::vector<GlobalPosition>;

public:
    // The default constructor
    FEIntegrationPointData() = delete;

    // The constructor
    template<class Geometry>
    FEIntegrationPointData(const Geometry& geometry,
                           const LocalPosition& ipLocal,
                           const LocalBasis& localBasis)
    : ipLocal_(ipLocal),
      ipGlobal_(geometry.global(ipLocal))
    {
        auto numLocalDofs = localBasis.size();

        // set the shape values
        shapeValues_.resize(numLocalDofs);
        localBasis.evaluateFunction(ipLocal, shapeValues_);

        // the local shape function gradients
        std::vector<JacobianType> shapeGrads(numLocalDofs);
        localBasis.evaluateJacobian(ipLocal, shapeGrads);

        // the global shape function gradients
        const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
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
    const LocalPosition& ipLocal() const
    { return ipLocal_; }

    //! The global position of the quadrature point
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

private:
    LocalPosition ipLocal_;
    GlobalPosition ipGlobal_;

    ShapeValues shapeValues_;
    ShapeGradients shapeGradients_;
};

/*!
 * \ingroup FEMDiscretization
 * \brief Integration point data related to a face of an element
 */
template<class GlobalPosition, class LocalBasis, class BoundaryFlag>
class FEFaceIntegrationPointData : public FEIntegrationPointData<GlobalPosition, LocalBasis>
{
    using ParentType = FEIntegrationPointData<GlobalPosition, LocalBasis>;
    using LocalPosition = typename LocalBasis::Traits::DomainType;
public:
    // The default constructor
    FEFaceIntegrationPointData() = delete;

    // The constructor
    template<class Geometry>
    FEFaceIntegrationPointData(const Geometry& geometry,
                               const LocalPosition& ipLocal,
                               const LocalBasis& localBasis,
                               const GlobalPosition& n,
                               const BoundaryFlag& bFlag)
    : ParentType(geometry, ipLocal, localBasis), normal_(n), boundaryFlag_(bFlag)
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
