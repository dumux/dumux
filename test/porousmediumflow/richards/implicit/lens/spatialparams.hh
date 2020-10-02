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
 * \ingroup RichardsTests
 * \brief Spatial parameters for the RichardsLensProblem.
 */

#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/fluidmatrix.hh>

#include <dumux/porousmediumflow/richards/model.hh>

namespace Dumux {

/*!
 * \ingroup RichardsTests
 * \brief The spatial parameters for the RichardsLensProblem.
 */
template<class GridGeometry, class Scalar>
class RichardsLensSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, RichardsLensSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RichardsLensSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum { dimWorld = GridView::dimensionworld };

public:
    using VanGenuchten = FluidMatrix::VanGenuchtenDefault<Scalar>;
    using PcSw = FluidMatrix::PcSw<FluidMatrix::VanGenuchtenDefault<Scalar>>;
    using KrSw = FluidMatrix::KrSw<FluidMatrix::VanGenuchtenDefault<Scalar>>;
    using FluidMatrixInteraction = FluidMatrix::FluidMatrixInteraction<PcSw, KrSw>;
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsLensSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {

        lensLowerLeft_ = {1.0, 2.0};
        lensUpperRight_ = {4.0, 3.0};

        // parameters for the Van Genuchten law
        // alpha and n
        typename FluidMatrixInteraction::PcSw::BasicParams lensParams(0.00045/*alpha*/, 7.3/*n*/);
        typename FluidMatrixInteraction::PcSw::BasicParams outerParams(0.0037/*alpha*/, 4.7/*n*/);

        // residual saturations
        typename FluidMatrixInteraction::PcSw::EffToAbsParams lensEffToAbsParams(0.18/*swr*/, 0.0/*snr*/);
        typename FluidMatrixInteraction::PcSw::EffToAbsParams outerEffToAbsParams(0.05/*swr*/, 0.0/*snr*/);

        auto vangGenuchtenLawLens = VanGenuchten(lensParams, lensEffToAbsParams);
        auto vangGenuchtenLawOuter = VanGenuchten(outerParams, outerEffToAbsParams);

        lensFluidMatrixInteraction_ = std::make_unique<FluidMatrixInteraction>(PcSw(vangGenuchtenLawLens), KrSw(vangGenuchtenLawLens));
        outerFluidMatrixInteraction_ = std::make_unique<FluidMatrixInteraction>(PcSw(vangGenuchtenLawOuter), KrSw(vangGenuchtenLawOuter));

        lensK_ = 1e-12;
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
    { return 0.4; }

    /*!
     * \brief Returns the parameters for the material law for the sub-control volume
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param element The current finite element
     * \param scv The sub-control volume
     * \param elemSol The current element solution
     */
    template<class ElementSolution>
    const FluidMatrixInteraction& fluidMatrixInteraction(const Element& element,
                                                         const SubControlVolume& scv,
                                                         const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        return fluidMatrixInteractionAtPos(globalPos);
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param globalPos The global coordinates for the given location
     */
    const FluidMatrixInteraction& fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return *lensFluidMatrixInteraction_;
        return *outerFluidMatrixInteraction_;
    }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] < lensLowerLeft_[i] - eps_ || globalPos[i] > lensUpperRight_[i] + eps_)
                return false;

        return true;
    }

    static constexpr Scalar eps_ = 1e-6;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;
    std::unique_ptr<FluidMatrixInteraction> lensFluidMatrixInteraction_;
    std::unique_ptr<FluidMatrixInteraction> outerFluidMatrixInteraction_;
};

} // end namespace Dumux

#endif
