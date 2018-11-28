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
 * \ingroup OnePTests
 * \brief The spatial params the incompressible test
 */
#ifndef DUMUX_ONEP_MILDANISOTROPY_SPATIAL_PARAMS_HH
#define DUMUX_ONEP_MILDANISOTROPY_SPATIAL_PARAMS_HH

#include <dune/common/fmatrix.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        incompressible 1p model
 */
template<class FVGridGeometry, class Scalar>
class MildAnisotropySpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             MildAnisotropySpatialParams<FVGridGeometry, Scalar>>
{
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           MildAnisotropySpatialParams<FVGridGeometry, Scalar>>;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = Tensor;

    MildAnisotropySpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        permeability_[0][0] = 1.5; permeability_[0][1] = 0.5;
        permeability_[1][0] = 0.5; permeability_[1][1] = 1.5;
    }

    //! Returns the permeability evaluated at a global position
    const PermeabilityType& permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    //! Returns the porosity evaluated at a global position
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    PermeabilityType permeability_;
};

} // end namespace Dumux

#endif
