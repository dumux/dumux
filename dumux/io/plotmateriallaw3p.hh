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
 * \brief Interface for plotting the two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_FLUID_MATRIX_LAW_HH
#define DUMUX_PLOT_FLUID_MATRIX_LAW_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/io/gnuplotinterface.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(MaterialLaw);
NEW_PROP_TAG(MaterialLawParams);
NEW_PROP_TAG(Scalar);
}

/*!
 *\brief Interface for plotting the three-phase fluid-matrix-interaction laws
 *
 * TODO: add theta head pressure plot (porosity and density is needed)
 */
template<class TypeTag>
class PlotMaterialLaw
{
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    //! Constructor
    PlotMaterialLaw(bool interaction = true)
    : numIntervals_(1000)
    {
        gnuplotpcgw_.setInteraction(interaction);
        gnuplotpcnw_.setInteraction(interaction);
        gnuplotpcgn_.setInteraction(interaction);
        gnuplotpcAlpha_.setInteraction(interaction);
        gnuplotdpcgw_dswe_.setInteraction(interaction);
        gnuplotdpcnw_dswe_.setInteraction(interaction);
        gnuplotdpcgn_dste_.setInteraction(interaction);
        gnuplotkr_.setInteraction(interaction);
        gnuplotkrn_.setInteraction(interaction);
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the water-gas interphase
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotpcgw(const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar pcMin = 0.0;
        Scalar pcMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName);

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pc[i] = MaterialLaw::pcgw(params, sw[i]);
            pcMin = std::min(pcMin, pc[i]);
            pcMax = std::max(pcMax, pc[i]);
        }

        gnuplotpcgw_.setXRange(lowerSat, upperSat);
        gnuplotpcgw_.setYRange(pcMin, pcMax);
        gnuplotpcgw_.setXlabel("wetting phase saturation [-]");
        gnuplotpcgw_.setYlabel("capillary pressure [Pa]");
        gnuplotpcgw_.addDataSetToPlot(sw, pc, plotName + "_pcgw-Sw");
        gnuplotpcgw_.plot("pcgw-Sw");
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the water-NAPL interface
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotpcnw(const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar pcMin = 0.0;
        Scalar pcMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName);

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pc[i] = MaterialLaw::pcnw(params, sw[i]);
            pcMin = std::min(pcMin, pc[i]);
            pcMax = std::max(pcMax, pc[i]);
        }

        gnuplotpcgw_.setXRange(lowerSat, upperSat);
        gnuplotpcgw_.setYRange(pcMin, pcMax);
        gnuplotpcgw_.setXlabel("wetting phase saturation [-]");
        gnuplotpcgw_.setYlabel("capillary pressure [Pa]");
        gnuplotpcgw_.addDataSetToPlot(sw, pc, plotName + "_pcnw-Sw");
        gnuplotpcgw_.plot("pcnw-Sw");
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the gas-NAPL interface
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotpcgn(const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "")
    {
        std::vector<Scalar> st(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar pcMin = 0.0;
        Scalar pcMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName);

        for (int i = 0; i <= numIntervals_; i++)
        {
            st[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pc[i] = MaterialLaw::pcgn(params, st[i]);
            pcMin = std::min(pcMin, pc[i]);
            pcMax = std::max(pcMax, pc[i]);
        }

        gnuplotpcgw_.setXRange(lowerSat, upperSat);
        gnuplotpcgw_.setYRange(pcMin, pcMax);
        gnuplotpcgw_.setXlabel("wetting phase saturation [-]");
        gnuplotpcgw_.setYlabel("capillary pressure [Pa]");
        gnuplotpcgw_.addDataSetToPlot(st, pc, plotName + "_pcgn-St");
        gnuplotpcgw_.plot("pcgn-St (St=Sw+Sn)");
    }


    /*!
     * \brief Plot the relative permeabilities
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotkr(const MaterialLawParams &params,
                Scalar lowerSat = 0.0,
                Scalar upperSat = 1.0,
                std::string plotName = "")
    {
        std::vector<Scalar> sw(numIntervals_ + 1);
        std::vector<Scalar> krw(numIntervals_ + 1);
        std::vector<Scalar> krg(numIntervals_ + 1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar krMin = 1e100;
        Scalar krMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName + "_kr");

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            krw[i] = MaterialLaw::krw(params, sw[i], 0.0);
            krg[i] = MaterialLaw::krg(params, sw[i], 0.0);
            krMin = std::min(krMin, std::min(krw[i], krg[i]));
            krMax = std::max(krMax, std::max(krw[i], krg[i]));
        }

        gnuplotkr_.setXRange(lowerSat, upperSat);
        gnuplotkr_.setYRange(krMin, krMax);
        gnuplotkr_.setXlabel("wetting phase saturation [-]");
        gnuplotkr_.setYlabel("relative permeability [-]");
        gnuplotkr_.addDataSetToPlot(sw, krw, plotName + "_krw");
        gnuplotkr_.addDataSetToPlot(sw, krg, plotName + "_krg");
        gnuplotkr_.plot("kr");
    }

    /*!
     * \brief Check the validity range for wetting saturation, to avoid an
     *        assert of the used material laws
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void checkEffectiveSaturation(const MaterialLawParams &params,
                                  Scalar lowerSat,
                                  Scalar upperSat,
                                  std::string plotName = "")
    {
        if (lowerSat < params.swr())
            Dune::dwarn << "warning: fluid-matrix law " << plotName << " can only be plotted for sw > swr" << std::endl;
        if (upperSat > (1.0 - params.snr()))
            Dune::dwarn << "warning: fluid-matrix law " << plotName << " can only be plotted for sw < 1.0 - snr" << std::endl;
    }

private:
    int numIntervals_;
    GnuplotInterface<Scalar> gnuplotpcgw_;
    GnuplotInterface<Scalar> gnuplotpcnw_;
    GnuplotInterface<Scalar> gnuplotpcgn_;
    GnuplotInterface<Scalar> gnuplotpcAlpha_;
    GnuplotInterface<Scalar> gnuplotdpcgw_dswe_;
    GnuplotInterface<Scalar> gnuplotdpcnw_dswe_;
    GnuplotInterface<Scalar> gnuplotdpcgn_dste_;
    GnuplotInterface<Scalar> gnuplotkr_;
    GnuplotInterface<Scalar> gnuplotkrn_;

};
} // end of namespace

#endif // DUMUX_PLOT_FLUID_MATRIX_LAW_HH
