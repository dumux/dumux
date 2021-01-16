// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Spatial parameters for the RichardsAnalyticalProblem.
 */

#ifndef DUMUX_RICHARDS_ANALYTICAL_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_ANALYTICAL_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

#include <dumux/porousmediumflow/richards/model.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsAnalyticalProblem.
 */
template<class GridGeometry, class Scalar>
class RichardsAnalyticalSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         RichardsAnalyticalSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       RichardsAnalyticalSpatialParams<GridGeometry, Scalar>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;

public:

    // export permeability type
    using PermeabilityType = Scalar;

    RichardsAnalyticalSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeability_ = 5e-12;

        typename PcKrSwCurve::BasicParams params(0/*pcEntry*/, 1e10/*pcMax*/);
        pcKrSwCurve_ = std::make_unique<PcKrSwCurve>(params);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     * \param globalPos The global position where we evaluate
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the porosity [] at a given location
     * \param globalPos The global position where we evaluate
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos A global coordinate vector
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return makeFluidMatrixInteraction(*pcKrSwCurve_);
    }

private:
    Scalar permeability_;
    std::unique_ptr<PcKrSwCurve> pcKrSwCurve_;
};

} // end namespace Dumux

#endif
