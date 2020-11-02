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
 * \ingroup ThreePTests
 * \brief Definition of the spatial parameters for the kuevette problem, which
 *        uses the three-phase fully implicit model.
 */

#ifndef DUMUX_INFILTRATION_THREEP_SPATIAL_PARAMS_HH
#define DUMUX_INFILTRATION_THREEP_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangenuchten.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/common/enumerate.hh>

namespace Dumux {
/*!
 * \ingroup ThreePTests
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class GridGeometry, class Scalar>
class InfiltrationThreePSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         InfiltrationThreePSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       InfiltrationThreePSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using ThreePhasePcKrSw = FluidMatrix::ParkerVanGenuchten3PDefault<Scalar>;

public:
    //! Export permeability type
    using PermeabilityType = Scalar;

    InfiltrationThreePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {
        // intrinsic permeabilities
        fineK_ = getParam<Scalar>("SpatialParams.permeability");
        coarseK_ = getParam<Scalar>("SpatialParams.permeability");

        // porosities
        porosity_ = getParam<Scalar>("SpatialParams.porosity");
        plotFluidMatrixInteractions_ =  getParam<bool>("Output.PlotFluidMatrixInteractions");
    }

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        GnuplotInterface<Scalar> gnuplot(plotFluidMatrixInteractions_);
        gnuplot.setOpenPlotWindow(plotFluidMatrixInteractions_);

        const Scalar sg = 0.2; // assume a fixed gas saturation
        auto swRange = linspace(0.2, 1.0, 1000);

        // assume fixed sn = 0.2 for pcgw curve
        auto pcgw = swRange;
        std::transform(swRange.begin(), swRange.end(), pcgw.begin(), [&](auto x){ return this->pcKrSwCurve_.pcgw(x, 0.2); });

        gnuplot.resetAll();
        gnuplot.setXlabel("Sw");
        gnuplot.setYlabel("pcgw");
        gnuplot.addDataSetToPlot(swRange, pcgw, "pcgw", "w l");
        gnuplot.plot("pcgw-sw");

        // plot kr
        swRange = linspace(0.2, 0.8, 1000);
        auto krw = swRange;
        auto krn = swRange;
        auto krg = swRange;
        for (const auto& [i, sw] : enumerate(swRange))
        {
            const Scalar sn = 1.0 - sg - sw;
            krw[i] = pcKrSwCurve_.krw(sw, sn);
            krn[i] = pcKrSwCurve_.krn(sw, sn);
            krg[i] = pcKrSwCurve_.krg(sw, sn);
        }

        gnuplot.resetAll();
        gnuplot.setXlabel("Sw");
        gnuplot.setYlabel("kr");
        gnuplot.addDataSetToPlot(swRange, krw, "krw", "w l");
        gnuplot.addDataSetToPlot(swRange, krn, "krn", "w l");
        gnuplot.addDataSetToPlot(swRange, krg, "krg", "w l");
        gnuplot.plot("kr-sw");
    }

      /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub-control volume
     * \param elemSol The element solution vector
     * \return The intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isFineMaterial_(scv.dofPosition()))
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
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    {
        return 70. - eps_ <= globalPos[0] && globalPos[0] <= 85. + eps_ &&
               7.0 - eps_ <= globalPos[1] && globalPos[1] <= 7.50 + eps_;
    }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar porosity_;
    const ThreePhasePcKrSw pcKrSwCurve_;

    bool plotFluidMatrixInteractions_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
