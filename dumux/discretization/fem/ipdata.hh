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
 * \ingroup FEMDiscretization
 * \brief A class to compute and store data on an
 *        integration point within a finite element.
 */
#ifndef DUMUX_DISCRETIZATION_FE_INTEGRATIONPOINT_DATA_HH
#define DUMUX_DISCRETIZATION_FE_INTEGRATIONPOINT_DATA_HH

#include <vector>

namespace Dumux {

/*!
 * \file
 * \ingroup FEMDiscretization
 * \brief Computes and stores data on an
 *        integration point within a finite element.
 * \tparam Geometry The geometry over which it is integrated
 * \tparam LocalBasis The local finite element basis
 */
template<class Geometry, class LocalBasis>
class FEIntegrationPointData
{
    using GlobalCoordinate = typename Geometry::GlobalCoordinate;
    using LocalCoordinate = typename LocalBasis::Traits::DomainType;
    using RangeType = typename LocalBasis::Traits::RangeType;
    using JacobianType = typename LocalBasis::Traits::JacobianType;

    using ShapeValues = std::vector<RangeType>;
    using ShapeGradients = std::vector<GlobalCoordinate>;

public:
    //! export the type used for global coordinates
    using GlobalPosition = GlobalCoordinate;

    /*!
     * \brief Constructor.
     * \param geometry the geometry of an element
     * \param ipLocal the local position of the integration point
     * \param localBasis the local finite element basis
     */
    FEIntegrationPointData(const Geometry& geometry,
                           const LocalCoordinate& ipLocal,
                           const LocalBasis& localBasis)
    : ipLocal_(ipLocal),
      ipGlobal_(geometry.global(ipLocal))
    {
        numLocalDofs_ = localBasis.size();

        // set the shape values
        shapeValues_.resize(numLocalDofs_);
        localBasis.evaluateFunction(ipLocal, shapeValues_);

        // the local shape function gradients
        std::vector<JacobianType> shapeGrads(numLocalDofs_);
        localBasis.evaluateJacobian(ipLocal, shapeGrads);

        // the global shape function gradients
        const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        shapeGradients_.resize(numLocalDofs_, GlobalPosition(0.0));
        for (unsigned int i = 0; i < numLocalDofs_; i++)
            jacInvT.umv(shapeGrads[i][0], shapeGradients_[i]);
    }

    //! Return the number of dofs for which data is stored
    std::size_t size() const
    { return numLocalDofs_; }

    //! The shape values at the quadrature point
    const ShapeValues& shapeValues() const
    { return shapeValues_; }

    //! The shape value of a local dof at the quadrature point
    const RangeType& shapeValue(unsigned int i) const
    { return shapeValues_[i]; }

    //! The shape value gradients at the quadrature point
    const ShapeGradients& shapeGradients() const
    { return shapeGradients_; }

    //! The shape value gradient of a local dof at the quadrature point
    //! \note the name of this is chosen such that it matches the one used in BoxFluxVariablesCache.
    const GlobalPosition& gradN(unsigned int i) const
    { return shapeGradients_[i]; }

    //! The local position of the quadrature point
    const LocalCoordinate& ipLocal() const
    { return ipLocal_; }

    //! The global position of the quadrature point
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

private:
    std::size_t numLocalDofs_;

    LocalCoordinate ipLocal_;
    GlobalPosition ipGlobal_;

    ShapeValues shapeValues_;
    ShapeGradients shapeGradients_;
};

} // end namespace Dumux

#endif
