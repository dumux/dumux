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
 * \brief TODO
 * \todo TODO Doc me.
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
        for (unsigned int i = 0; i < numLocalDofs; i++)
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
    //! \note the name of this interface is chosen such that it matches the
    //!       one used in BoxFluxVariablesCache.
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

} // end namespace

#endif
