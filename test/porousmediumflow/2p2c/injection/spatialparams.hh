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
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component fully implicit model.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component fully implicit model.
 */
template<class GridGeometry, class Scalar>
class InjectionSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         InjectionSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, InjectionSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;

    InjectionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , finePcKrSwCurve_("SpatialParams.FineMaterial")
    , coarsePcKrSwCurve_("SpatialParams.CoarseMaterial")
    {
        layerBottom_ = 22.5;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        return coarsePorosity_;
    }


    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return makeFluidMatrixInteraction(finePcKrSwCurve_);
        return makeFluidMatrixInteraction(coarsePcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The position of the center of the element
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > layerBottom_; }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    const PcKrSwCurve finePcKrSwCurve_;
    const PcKrSwCurve coarsePcKrSwCurve_;
};

} // end namespace Dumux

#endif
