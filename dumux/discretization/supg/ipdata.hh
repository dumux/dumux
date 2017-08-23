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
 * \brief Base class for a sub control volume
 */
#ifndef DUMUX_DISCRETIZATION_FEM_INTEGRATION_POINT_DATA_HH
#define DUMUX_DISCRETIZATION_FEM_INTEGRATION_POINT_DATA_HH

#include <dune/common/fvector.hh>

namespace Dumux
{

template<class TypeTag>
class FemIntegrationPointData
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FEBasis = typename GET_PROP_TYPE(TypeTag, FeBasis);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using Element = typename GridView::template Codim<0>::Entity;
    using Geometry = typename Element::Geometry;
    using LocalView = typename FEBasis::LocalView;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using LocalPosition = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using ShapeValue = Dune::FieldVector<Scalar, 1>;
    using ShapeJacobian = Dune::FieldMatrix<Scalar, 1, dim>;
    using ShapeValues = std::vector<ShapeValue>;
    using RefShapeGradients = std::vector<ShapeJacobian>;
    using ShapeGradients = std::vector<GlobalPosition>;

public:
    // The default constructor
    FemIntegrationPointData() = delete;

    // The constructor
    template<class LocalBasis>
    FemIntegrationPointData(const Geometry& geometry,
                            const LocalPosition& ipLocal,
                            const LocalBasis& localBasis)
    : ipLocal_(ipLocal),
      ipGlobal_(geometry.global(ipLocal))
    {
        auto numLocalDofs = localBasis.size();

        // set the shape values
        shapeValues_.resize(numLocalDofs);
        localBasis.evaluateFunction(ipLocal, shapeValues_);

        // the shape gradients on the reference element
        RefShapeGradients shapeGradRef(numLocalDofs);
        localBasis.evaluateJacobian(ipLocal, shapeGradRef);

        // the global shape function gradients
        const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        shapeGradients_.resize(numLocalDofs, GlobalPosition(0.0));
        for (size_t i = 0; i < numLocalDofs; i++)
            jacInvT.umv(shapeGradRef[i][0],shapeGradients_[i]);
    }

    //! The shape values at the quadrature point
    const ShapeValues& shapeValues() const
    { return shapeValues_; }

    //! The shape value of local dof at the quadrature point
    const ShapeValue& shapeValues(int i) const
    { return shapeValues_[i]; }

    //! The shape value gradients at the quadrature point
    const ShapeGradients& shapeGradients() const
    { return shapeGradients_; }

    //! The shape value gradient of local dof at the quadrature point
    const GlobalPosition& shapeGradients(int i) const
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
