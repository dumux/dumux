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
 * \ingroup RichardsNCTests
 * \brief The spatial parameters for the RichardsWellTracerProblem.
 */

#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

#include <dumux/porousmediumflow/richards/model.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsWellTracerProblem.
 */
template<class GridGeometry, class Scalar>
class RichardsWellTracerSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         RichardsWellTracerSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       RichardsWellTracerSpatialParams<GridGeometry, Scalar>>;

    enum { dimWorld=GridView::dimensionworld };
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenNoReg<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsWellTracerSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , lensPcKrSwCurve_("SpatialParams.Lens")
    , outerPcKrSwCurve_("SpatialParams.OuterDomain")
    {
        lensLowerLeft_ = getParam<GlobalPosition>("Problem.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("Problem.LensUpperRight");

        lensK_ = 1e-14;
        outerK_ = 5e-12;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return 0.2;
        return 0.4;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return makeFluidMatrixInteraction(lensPcKrSwCurve_);
        else
            return makeFluidMatrixInteraction(outerPcKrSwCurve_);
    }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] < lensLowerLeft_[i] - eps_ || globalPos[i] > lensUpperRight_[i] + eps_)
                return false;

        return true;
    }

    static constexpr Scalar eps_ = 1.5e-7;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;

    const PcKrSwCurve lensPcKrSwCurve_;
    const PcKrSwCurve outerPcKrSwCurve_;
};

} // end namespace Dumux

#endif
