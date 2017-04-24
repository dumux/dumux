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
 * \brief Definition of the spatial parameters for the effective diffusivity tests.
 */
#ifndef DUMUX_FLUIDMATRIXINTERACTION_TEST_SPATIAL_PARAMS_TWOP_HH
#define DUMUX_FLUIDMATRIXINTERACTION_TEST_SPATIAL_PARAMS_TWOP_HH

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>

#include "../fluidmatrixinteractionsspatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class ThermalConductivityTestSpatialParamsTwoP
 : public FluidMatrixInteractionTestSpatialParams<TypeTag>
{
    using ParentType = FluidMatrixInteractionTestSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    ThermalConductivityTestSpatialParamsTwoP(const Problem& problem, const GridView &gridView)
    : ParentType(problem, gridView) {}
    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        GnuplotInterface<Scalar> gnuplot;
        gnuplot.setOpenPlotWindow(GET_PARAM_FROM_GROUP(TypeTag, bool, Output, OpenPlotWindow));
        PlotThermalConductivityModel<TypeTag> plotThermalConductivityModel_(293.15, 1e5);
        std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, File);
        plotThermalConductivityModel_.addlambdaeffcurve(gnuplot, this->porosity_, this->rhoSolid_, this->lambdaSolid_,
                                                        0.0, 1.0, fileName + ".dat");
        gnuplot.plot("lambda_eff");
    }
};

} // end namespace Dumux

#endif
