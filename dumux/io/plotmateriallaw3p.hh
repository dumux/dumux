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
 * \brief Interface for plotting the three-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_FLUID_MATRIX_LAW_HH
#define DUMUX_PLOT_FLUID_MATRIX_LAW_HH

#include <dune/common/deprecated.hh>

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
        gnuplotpc_.setInteraction(interaction);
        gnuplotpcAlpha_.setInteraction(interaction);
        gnuplotkr_.setInteraction(interaction);
        gnuplotkrn_.setInteraction(interaction);
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for all  phases
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotpc(const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "")
    {
        plotpcgw(params, lowerSat, upperSat, plotName);
        plotpcnw(params, lowerSat, upperSat, plotName);
        plotpcgn(params, lowerSat, upperSat, plotName);
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
        std::vector<Scalar> sw;
        std::vector<Scalar> pc;
        Scalar satInterval = upperSat - lowerSat;
        Scalar pcMin = 0.0;
        Scalar pcMax = -1e100;

        Scalar swTemp, pcTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pcTemp = MaterialLaw::pcgw(params, swTemp);
            if (checkValues(swTemp, pcTemp))
            {
                sw.push_back(swTemp);
                pc.push_back(pcTemp);
                pcMin = std::min(pcMin, pcTemp);
                pcMax = std::max(pcMax, pcTemp);
            }
        }

        gnuplotpc_.setXRange(lowerSat, upperSat);
        gnuplotpc_.setYRange(pcMin, pcMax);
        gnuplotpc_.setXlabel("wetting phase saturation [-]");
        gnuplotpc_.setYlabel("capillary pressure [Pa]");
        gnuplotpc_.addDataSetToPlot(sw, pc, plotName + "_pcgw-Sw");
        gnuplotpc_.plot("pcgw-Sw");
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
        std::vector<Scalar> sw;
        std::vector<Scalar> pc;
        Scalar satInterval = upperSat - lowerSat;
        Scalar pcMin = 0.0;
        Scalar pcMax = -1e100;

        Scalar swTemp, pcTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pcTemp = MaterialLaw::pcnw(params, swTemp);
            if (checkValues(swTemp, pcTemp))
            {
                sw.push_back(swTemp);
                pc.push_back(pcTemp);
                pcMin = std::min(pcMin, pcTemp);
                pcMax = std::max(pcMax, pcTemp);
            }
        }

        gnuplotpc_.setXRange(lowerSat, upperSat);
        gnuplotpc_.setYRange(pcMin, pcMax);
        gnuplotpc_.setXlabel("wetting phase saturation [-]");
        gnuplotpc_.setYlabel("capillary pressure [Pa]");
        gnuplotpc_.addDataSetToPlot(sw, pc, plotName + "_pcnw-Sw");
        gnuplotpc_.plot("pcnw-Sw");
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
        std::vector<Scalar> st;
        std::vector<Scalar> pc;
        Scalar satInterval = upperSat - lowerSat;
        Scalar pcMin = 0.0;
        Scalar pcMax = -1e100;

        Scalar stTemp, pcTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            stTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pcTemp = MaterialLaw::pcgn(params, stTemp);
            if (checkValues(stTemp, pcTemp))
            {
                st.push_back(stTemp);
                pc.push_back(pcTemp);
                pcMin = std::min(pcMin, pcTemp);
                pcMax = std::max(pcMax, pcTemp);
            }
        }

        gnuplotpc_.setXRange(lowerSat, upperSat);
        gnuplotpc_.setYRange(pcMin, pcMax);
        gnuplotpc_.setXlabel("wetting phase saturation [-]");
        gnuplotpc_.setYlabel("capillary pressure [Pa]");
        gnuplotpc_.addDataSetToPlot(st, pc, plotName + "_pcgn-St");
        gnuplotpc_.plot("pcgn-St (St=Sw+Sn)");
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
        std::vector<Scalar> krn(numIntervals_ + 1);
        std::vector<Scalar> krg(numIntervals_ + 1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar krMin = 1e100;
        Scalar krMax = -1e100;

        Scalar swTemp, krwTemp, krnTemp, krgTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            krwTemp = MaterialLaw::krw(params, swTemp);
            krnTemp = MaterialLaw::krn(params, swTemp);
            krgTemp = MaterialLaw::krg(params, swTemp);
            if (checkValues(swTemp, krwTemp)
                && checkValues(swTemp, krnTemp)
                && checkValues(swTemp, krgTemp))
            {
                sw.push_back(swTemp);
                krw.push_back(krwTemp);
                krn.push_back(krnTemp);
                krn.push_back(krgTemp);
                krMin = std::min({krMin, krwTemp, krnTemp, krgTemp});
                krMax = std::max({krMax, krwTemp, krnTemp, krgTemp});
            }
        }

        gnuplotkr_.setXRange(lowerSat, upperSat);
        gnuplotkr_.setYRange(krMin, krMax);
        gnuplotkr_.setXlabel("wetting phase saturation [-]");
        gnuplotkr_.setYlabel("relative permeability [-]");
        gnuplotkr_.addDataSetToPlot(sw, krw, plotName + "_krw");
        gnuplotkr_.addDataSetToPlot(sw, krn, plotName + "_krn");
        gnuplotkr_.addDataSetToPlot(sw, krg, plotName + "_krg");
        gnuplotkr_.plot("kr");
    }

    /*!
     * \brief Plot the transition (2P/3P) function
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotPcAlpha(const MaterialLawParams &params,
                Scalar lowerSat = 0.0,
                Scalar upperSat = 1.0,
                std::string plotName = "")
    {
        std::vector<Scalar> sn(numIntervals_ + 1);
        std::vector<Scalar> alpha(numIntervals_ + 1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar alphaMin = -2;
        Scalar alphaMax = 2;

        Scalar snTemp, alphaTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            snTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            alphaTemp = MaterialLaw::pcAlpha(params, snTemp);
            if (checkValues(snTemp, alphaTemp))
            {
                sn.push_back(snTemp);
                alpha.push_back(alphaTemp);
                alphaMin = std::min(alphaMin, alphaTemp);
                alphaMax = std::max(alphaMax, alphaTemp);
            }
        }

        gnuplotpcAlpha_.setXRange(lowerSat, upperSat);
        gnuplotpcAlpha_.setYRange(alphaMin, alphaMax);
        gnuplotpcAlpha_.setXlabel("non-wetting phase saturation [-]");
        gnuplotpcAlpha_.setYlabel("transition function [-]");
        gnuplotpcAlpha_.addDataSetToPlot(sn, alpha, plotName + "_alpha");
        gnuplotpcAlpha_.plot("alpha");
    }

private:
    /*!
     * \brief Check the validity range for wetting saturation, to avoid an
     *        assert of the used material laws
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    DUNE_DEPRECATED_MSG("checkEffectiveSaturation is deprecated.")
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

    /*!
     * \brief Check the values for occurrences of nan and inf
     *
     * \param value1 A data point value
     * \param value2 An other data point value
     */
    bool checkValues(Scalar value1, Scalar value2)
    {
        return !std::isnan(value1) && !std::isinf(value1)
               && !std::isnan(value2) && !std::isinf(value2);
    }

    int numIntervals_;
    GnuplotInterface<Scalar> gnuplotpc_;
    GnuplotInterface<Scalar> gnuplotpcAlpha_;
    GnuplotInterface<Scalar> gnuplotkr_;
    GnuplotInterface<Scalar> gnuplotkrn_;

};
} // end of namespace

#endif // DUMUX_PLOT_FLUID_MATRIX_LAW_HH
