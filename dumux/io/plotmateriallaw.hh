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
 *\brief Interface for plotting the two-phase fluid-matrix-interaction laws
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
        gnuplotpcsw_.setInteraction(interaction);
        gnuplotswpc_.setInteraction(interaction);
        gnuplotdpcdsw_.setInteraction(interaction);
        gnuplotdswdpc_.setInteraction(interaction);
        gnuplotkr_.setInteraction(interaction);
        gnuplotkrdsw_.setInteraction(interaction);
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void plotpcsw(const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string curveTitle = "")
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
            pcTemp = MaterialLaw::pc(params, swTemp);
            if (checkValues_(swTemp, pcTemp))
            {
                sw.push_back(swTemp);
                pc.push_back(pcTemp);
                pcMin = std::min(pcMin, pcTemp);
                pcMax = std::max(pcMax, pcTemp);
            }
        }

        // use log scale for very high capillary pressures
        if (pcMax / std::max(pcMin, 1.0) > 1e6)
        {
            gnuplotpcsw_.setOption("set log y");
            pcMin = std::max(pcMin, 1.0);
        }

        gnuplotpcsw_.setXRange(lowerSat, upperSat);
        gnuplotpcsw_.setYRange(pcMin, pcMax);
        gnuplotpcsw_.setXlabel("wetting phase saturation [-]");
        gnuplotpcsw_.setYlabel("capillary pressure [Pa]");
        gnuplotpcsw_.addDataSetToPlot(sw, pc, curveTitle + "_pc-Sw");
        gnuplotpcsw_.plot("pc-Sw");
    }

    /*!
     * \brief Plot the saturation-capillary pressure curve
     *
     * \param params The material law parameters
     * \param lowerpc Minimum x-value
     * \param upperpc Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void plotswpc(const MaterialLawParams &params,
                  Scalar lowerpc = 0.0,
                  Scalar upperpc = 5000.0,
                  std::string curveTitle = "")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> pc;
        Scalar pcInterval = upperpc - lowerpc;
        Scalar swMin = 1e100;
        Scalar swMax = -1e100;

        Scalar pcTemp, swTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            pcTemp = lowerpc + pcInterval * Scalar(i) / Scalar(numIntervals_);
            swTemp = MaterialLaw::sw(params, pcTemp);
            if (checkValues_(pcTemp, swTemp))
            {
                pc.push_back(pcTemp);
                sw.push_back(swTemp);
                swMin = std::min(swMin, swTemp);
                swMax = std::max(swMax, swTemp);
            }
        }

        gnuplotswpc_.setXRange(lowerpc, upperpc);
        gnuplotswpc_.setYRange(swMin, swMax);
        gnuplotswpc_.setXlabel("capillary pressure [Pa]");
        gnuplotswpc_.setYlabel("wetting phase saturation [-]");
        gnuplotswpc_.addDataSetToPlot(pc, sw, curveTitle + "_Sw-pc");
        gnuplotswpc_.plot("sw-pc");
    }

    /*!
     * \brief Plot the gradient of the capillary pressure-saturation curve
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void plotdpcdsw(const MaterialLawParams &params,
                    Scalar lowerSat = 0.0,
                    Scalar upperSat = 1.0,
                    std::string curveTitle = "")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> dpcdsw;
        Scalar satInterval = upperSat - lowerSat;
        Scalar dpcdswMin = 1e100;
        Scalar dpcdswMax = -1e100;

        Scalar swTemp, dpcdswTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            dpcdswTemp = MaterialLaw::dpc_dsw(params, swTemp);
            if (checkValues_(swTemp, dpcdsw))
            {
                sw.push_back(swTemp);
                dpcdsw.push_back(dpcdswTemp);
                dpcdswMin = std::min(dpcdswMin, dpcdswTemp);
                dpcdswMax = std::max(dpcdswMax, dpcdswTemp);
            }
        }

        gnuplotdpcdsw_.setXRange(lowerSat, upperSat);
        gnuplotdpcdsw_.setYRange(dpcdswMin, dpcdswMax);
        gnuplotdpcdsw_.setXlabel("wetting phase saturation [-]");
        gnuplotdpcdsw_.setYlabel("gradient of the pc-Sw curve [Pa]");
        gnuplotdpcdsw_.addDataSetToPlot(sw, dpcdsw, curveTitle + "_dpcdSw-Sw");
        gnuplotdpcdsw_.plot("dpcdsw");
    }

    /*!
     * \brief Plot the gradient of the saturation-capillary pressure curve
     *
     * \param params The material law parameters
     * \param lowerpc Minimum x-value
     * \param upperpc Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void plotdswdpc(const MaterialLawParams &params,
                    Scalar lowerpc = 0.0,
                    Scalar upperpc = 5000.0,
                    std::string curveTitle = "")
    {
        std::vector<Scalar> pc;
        std::vector<Scalar> dswdpc;
        Scalar pcInterval = upperpc - lowerpc;
        Scalar dswdpcMin = 1e100;
        Scalar dswdpcMax = -1e100;

        Scalar dswdpcTemp, pcTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            pcTemp = lowerpc + pcInterval * Scalar(i) / Scalar(numIntervals_);
            dswdpcTemp = MaterialLaw::dsw_dpc(params, pcTemp);
            if (checkValues_(pcTemp, dswdpcTemp))
            {
                pc.push_back(pcTemp);
                dswdpc.push_back(dswdpcTemp);
                dswdpcMin = std::min(dswdpcMin, dswdpcTemp);
                dswdpcMax = std::max(dswdpcMax, dswdpcTemp);
            }
        }

        gnuplotdswdpc_.setXRange(lowerpc, upperpc);
        gnuplotdswdpc_.setYRange(dswdpcMin, dswdpcMax);
        gnuplotdswdpc_.setXlabel("capillary pressure [Pa]");
        gnuplotdswdpc_.setYlabel("gradient of the Sw-pc curve [1/Pa]");
        gnuplotdswdpc_.addDataSetToPlot(pc, dswdpc, curveTitle + "_dSwdpc-pc");
        gnuplotdswdpc_.plot("dswdpc");
    }

    /*!
     * \brief Plot the relative permeabilities
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void plotkr(const MaterialLawParams &params,
                Scalar lowerSat = 0.0,
                Scalar upperSat = 1.0,
                std::string curveTitle = "")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> krw;
        std::vector<Scalar> krn;
        Scalar satInterval = upperSat - lowerSat;
        Scalar krMin = 1e100;
        Scalar krMax = -1e100;

        Scalar swTemp, krwTemp, krnTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            krwTemp = MaterialLaw::krw(params, swTemp);
            krnTemp = MaterialLaw::krn(params, swTemp);
            if (checkValues_(swTemp, krwTemp) && checkValues_(swTemp, krnTemp))
            {
                sw.push_back(swTemp);
                krw.push_back(krwTemp);
                krn.push_back(krnTemp);
                krMin = std::min({krMin, krwTemp, krnTemp});
                krMax = std::max({krMax, krwTemp, krnTemp});
            }
        }

        gnuplotkr_.setXRange(lowerSat, upperSat);
        gnuplotkr_.setYRange(krMin, krMax);
        gnuplotkr_.setXlabel("wetting phase saturation [-]");
        gnuplotkr_.setYlabel("relative permeability [-]");
        gnuplotkr_.addDataSetToPlot(sw, krw, curveTitle + "_krw");
        gnuplotkr_.addDataSetToPlot(sw, krn, curveTitle + "_krn");
        gnuplotkr_.plot("kr");
    }

    /*!
     * \brief Plot the gradient of the relative permeabilities
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void plotdkrdsw(const MaterialLawParams &params,
                    Scalar lowerSat = 0.0,
                    Scalar upperSat = 1.0,
                    std::string curveTitle = "")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> dkrw_dsw;
        std::vector<Scalar> dkrn_dsw;
        Scalar satInterval = upperSat - lowerSat;
        Scalar dkrdswMin = 1e100;
        Scalar dkrdswMax = -1e100;

        Scalar swTemp, dkrwdswTemp, dkrndswTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            dkrwdswTemp = MaterialLaw::dkrw_dsw(params, swTemp);
            dkrndswTemp = MaterialLaw::dkrn_dsw(params, swTemp);
            if (checkValues_(swTemp, dkrwdswTemp) && checkValues_(swTemp, dkrndswTemp))
            {
                sw.push_back(swTemp);
                dkrw_dsw.push_back(dkrwdswTemp);
                dkrn_dsw.push_back(dkrndswTemp);
                dkrdswMin = std::min({dkrdswMin, dkrwdswTemp, dkrndswTemp});
                dkrdswMax = std::max({dkrdswMax, dkrwdswTemp, dkrndswTemp});
            }
        }

        gnuplotkrdsw_.setXRange(lowerSat, upperSat);
        gnuplotkrdsw_.setYRange(dkrdswMin, dkrdswMax);
        gnuplotkrdsw_.setXlabel("wetting phase saturation [-]");
        gnuplotkrdsw_.setYlabel("gradient of the kr-Sw function [-]");
        gnuplotkrdsw_.addDataSetToPlot(sw, dkrw_dsw, curveTitle + "_dkrw_dsw");
        gnuplotkrdsw_.addDataSetToPlot(sw, dkrn_dsw, curveTitle + "_dkrn_dsw");
        gnuplotkrdsw_.plot("dkrndsw");
    }

private:
    /*!
     * \brief Check the values for occurrences of nan and inf
     *
     * \param value1 A data point value
     * \param value2 An other data point value
     */
    bool checkValues_(Scalar value1, Scalar value2)
    {
        return !std::isnan(value1) && !std::isinf(value1)
               && !std::isnan(value2) && !std::isinf(value2);
    }

    int numIntervals_;
    GnuplotInterface<Scalar> gnuplotpcsw_;
    GnuplotInterface<Scalar> gnuplotswpc_;
    GnuplotInterface<Scalar> gnuplotdpcdsw_;
    GnuplotInterface<Scalar> gnuplotdswdpc_;
    GnuplotInterface<Scalar> gnuplotkr_;
    GnuplotInterface<Scalar> gnuplotkrdsw_;
};
} // end of namespace

#endif // DUMUX_PLOT_FLUID_MATRIX_LAW_HH
