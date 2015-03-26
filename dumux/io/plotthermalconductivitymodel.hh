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
 *
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH
#define DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(ThermalConductivityModel);
NEW_PROP_TAG(CompositionalFluidState);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(Scalar);
}

/*!
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
template<class TypeTag>
class PlotThermalConductivityModel
{
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

public:
    //! Constructor, initializes the fluid system
    PlotThermalConductivityModel()
    : numIntervals_(1000)
    {
        FluidState fluidstate;
        lambdaW_ = FluidSystem::template thermalConductivity<FluidState>(fluidstate, wPhaseIdx);
        lambdaN_ = FluidSystem::template thermalConductivity<FluidState>(fluidstate, nPhaseIdx);
    }

    /*!
     * \brief Plot the effective thermal conductivity-saturation curve
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotlambdaeff(Scalar porosity, Scalar rhoSolid, Scalar lambdaSolid,
                       Scalar lowerSat = 0.0,
                       Scalar upperSat = 1.0,
                       std::string plotName = "lambda_eff",
                       bool interaction = true)
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> lambda(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar lambdaMin = 0.0;
        Scalar lambdaMax = -1e100;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            lambda[i] = ThermalConductivityModel::effectiveThermalConductivity(sw[i], lambdaW_,
                                                                               lambdaN_, lambdaSolid,
                                                                               porosity, rhoSolid);
            lambdaMin = std::min(lambdaMin, lambda[i]);
            lambdaMax = std::max(lambdaMax, lambda[i]);
        }

        unsigned int windowNumber = 8;
        gnuplot_.setXRange(lowerSat, upperSat,windowNumber);
        gnuplot_.setYRange(lambdaMin, lambdaMax, windowNumber);
        gnuplot_.setXlabel("wetting phase saturation [-]", windowNumber);
        gnuplot_.setYlabel("effective thermal conductivity [W/(m K)]", windowNumber);
        gnuplot_.addDataSetToPlot(sw, lambda, plotName, windowNumber);
        gnuplot_.plot("lambdaeff", windowNumber, interaction);
    }

private:
    GnuplotInterface<Scalar> gnuplot_;
    int numIntervals_;
    Scalar lambdaN_;
    Scalar lambdaW_;
};
} // end of namespace

#endif // DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH
