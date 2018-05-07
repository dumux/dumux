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
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
#ifndef DUMUX_FRACFLOW_CONSTVEL_SPATIAL_PARAMS_HH
#define DUMUX_FRACFLOW_CONSTVEL_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedtestlaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/testlaw.hh>
#include <dumux/io/plotmateriallaw.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief The spatial parameters for the ConstVelProblem which uses the
 *        two-phase fully implicit model
 */
template<class FVGridGeometry, class Scalar>
class ConstVelSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         ConstVelSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = ConstVelSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
public:
    // export permeability type
    using PermeabilityType = Scalar;

    // export material law type
    using MaterialLaw = EffToAbsLaw<Testlaw<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    ConstVelSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // residual saturations
        materialParams_.setSwr(0.0);
        materialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        // alpha and n
        materialParams_.setPe(5171.068);
        materialParams_.setLambda(4.0);

        // regularization
        // materialParams_.setThresholdSw(getParam<Scalar>("Problem.RegularizationSw"));

        K_ = 9.86923* 1.0e-13; // 1 Darcy = 1000 milliDarcy

        if (getParam<bool>("Problem.PlotMaterialLaw"))
        {
            static const auto plotRange = getParam<std::vector<Scalar>>("Problem.PlotRange");
            static const std::string rangeString = "[" + std::to_string(plotRange[0]) + ":" + std::to_string(plotRange[1]) + "]";
            PlotMaterialLaw<Scalar, MaterialLaw> plotMaterialLaw;

            GnuplotInterface<Scalar> gnuplot(true);
            plotMaterialLaw.addDswcurve( gnuplot, materialParams_, plotRange[0], plotRange[1]);
            gnuplot.setOption("set xrange " + rangeString);
            gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
            gnuplot.plot("D-Sw");

            GnuplotInterface<Scalar> gnuplot2(true);
            plotMaterialLaw.addpcswcurve( gnuplot2, materialParams_, plotRange[0], plotRange[1]);
            gnuplot2.setOption("set xrange " + rangeString);
            gnuplot2.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
            gnuplot2.plot("pc-Sw");

            GnuplotInterface<Scalar> gnuplot3(true);
            plotMaterialLaw.addkrcurves( gnuplot3, materialParams_, plotRange[0], plotRange[1]);
            gnuplot3.setOption("set xrange " + rangeString);
            gnuplot3.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
            gnuplot3.plot("krw-Sw");

            GnuplotInterface<Scalar> gnuplot4(true);
            plotMaterialLaw.addkrcurves( gnuplot4, materialParams_, plotRange[0], plotRange[1]);
            gnuplot4.setOption("set xrange " + rangeString);
            gnuplot4.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
            gnuplot4.plot("kr-Sw");

            GnuplotInterface<Scalar> gnuplot5(true);
            plotMaterialLaw.adddpcdswcurve( gnuplot5, materialParams_, plotRange[0], plotRange[1]);
            gnuplot5.setOption("set xrange " + rangeString);
            gnuplot5.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
            gnuplot5.plot("dpc-dSw");

            exit(0);
        }
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    const PermeabilityType& permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return K_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx; // water phase is wetting
    }

private:
    Scalar K_;
    MaterialLawParams materialParams_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
