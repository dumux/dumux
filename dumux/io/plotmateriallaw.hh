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
    PlotMaterialLaw()
    : numIntervals_(1000)
    { }

    /*!
     * \brief Plot the capillary pressure-saturation curve
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotpcsw(const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "pc-Sw",
                  bool interaction = true)
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
            pc[i] = MaterialLaw::pc(params, sw[i]);
            pcMin = std::min(pcMin, pc[i]);
            pcMax = std::max(pcMax, pc[i]);
        }

        unsigned int windowNumber = 1;
        gnuplot_.setXRange(lowerSat, upperSat, windowNumber);
        gnuplot_.setYRange(pcMin, pcMax, windowNumber);
        gnuplot_.setXlabel("wetting phase saturation [-]", windowNumber);
        gnuplot_.setYlabel("capillary pressure [Pa]", windowNumber);
        gnuplot_.addDataSetToPlot(sw, pc, plotName, windowNumber);
        gnuplot_.plot("pc-Sw", windowNumber, interaction);
    }

    /*!
     * \brief Plot the saturation-capillary pressure curve
     *
     * \param params The material law parameters
     * \param lowerpc Minimum x-value
     * \param upperpc Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotswpc(const MaterialLawParams &params,
                  Scalar lowerpc = 0.0,
                  Scalar upperpc = 5000.0,
                  std::string plotName = "Sw-pc",
                  bool interaction = true)
    {
        std::vector<Scalar> sat(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        Scalar pcInterval = upperpc - lowerpc;
        Scalar swMin = 1e100;
        Scalar swMax = -1e100;

        for (int i = 0; i <= numIntervals_; i++)
        {
            pc[i] = lowerpc + pcInterval * Scalar(i) / Scalar(numIntervals_);
            sat[i] = MaterialLaw::sw(params, pc[i]);
            swMin = std::min(swMin, sat[i]);
            swMax = std::max(swMax, sat[i]);
        }

        unsigned int windowNumber = 2;
        gnuplot_.setXRange(lowerpc, upperpc, windowNumber);
        gnuplot_.setYRange(swMin, swMax, windowNumber);
        gnuplot_.setXlabel("capillary pressure [Pa]", windowNumber);
        gnuplot_.setYlabel("wetting phase saturation [-]", windowNumber);
        gnuplot_.addDataSetToPlot(pc, sat, plotName, windowNumber);
        gnuplot_.plot("sw-pc", windowNumber, interaction);
    }

    /*!
     * \brief Plot the gradient of the capillary pressure-saturation curve
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotdpcdsw(const MaterialLawParams &params,
                    Scalar lowerSat = 0.0,
                    Scalar upperSat = 1.0,
                    std::string plotName = "dpcdSw-Sw",
                    bool interaction = true)
    {
        std::vector<Scalar> sat(numIntervals_ + 1);
        std::vector<Scalar> dpcdsw(numIntervals_ + 1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar dpcdswMin = 1e100;
        Scalar dpcdswMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName);

        for (int i = 0; i <= numIntervals_; i++)
        {
            sat[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            dpcdsw[i] = MaterialLaw::dpc_dsw(params, sat[i]);
            dpcdswMin = std::min(dpcdswMin, dpcdsw[i]);
            dpcdswMax = std::max(dpcdswMax, dpcdsw[i]);
        }

        unsigned int windowNumber = 3;
        gnuplot_.setXRange(lowerSat, upperSat, windowNumber);
        gnuplot_.setYRange(dpcdswMin, dpcdswMax, windowNumber);
        gnuplot_.setXlabel("wetting phase saturation [-]", windowNumber);
        gnuplot_.setYlabel("gradient of the pc-Sw curve [Pa]", windowNumber);
        gnuplot_.addDataSetToPlot(sat, dpcdsw, plotName, windowNumber);
        gnuplot_.plot("dpcdsw", windowNumber, interaction);
    }

    /*!
     * \brief Plot the gradient of the saturation-capillary pressure curve
     *
     * \param params The material law parameters
     * \param lowerpc Minimum x-value
     * \param upperpc Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotdswdpc(const MaterialLawParams &params,
                    Scalar lowerpc = 0.0,
                    Scalar upperpc = 5000.0,
                    std::string plotName = "dSwdpc-pc",
                    bool interaction = true)
    {
        std::vector<Scalar> dswdpc(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        Scalar pcInterval = upperpc - lowerpc;
        Scalar dswdpcMin = 1e100;
        Scalar dswdpcMax = -1e100;

        for (int i = 0; i <= numIntervals_; i++)
        {
            pc[i] = lowerpc + pcInterval * Scalar(i) / Scalar(numIntervals_);
            dswdpc[i] = MaterialLaw::dsw_dpc(params, pc[i]);
            dswdpcMin = std::min(dswdpcMin, dswdpc[i]);
            dswdpcMax = std::max(dswdpcMax, dswdpc[i]);
        }

        unsigned int windowNumber = 4;
        gnuplot_.setXRange(lowerpc, upperpc, windowNumber);
        gnuplot_.setYRange(dswdpcMin, dswdpcMax, windowNumber);
        gnuplot_.setXlabel("capillary pressure [Pa]", windowNumber);
        gnuplot_.setYlabel("gradient of the Sw-pc curve [1/Pa]", windowNumber);
        gnuplot_.addDataSetToPlot(pc, dswdpc, plotName, windowNumber);
        gnuplot_.plot("dswdpc", windowNumber, interaction);
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
                std::string plotName = "",
                bool interaction = true)
    {
        std::vector<Scalar> sw(numIntervals_ + 1);
        std::vector<Scalar> krw(numIntervals_ + 1);
        std::vector<Scalar> krn(numIntervals_ + 1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar krMin = 1e100;
        Scalar krMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName + "_kr");

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            krw[i] = MaterialLaw::krw(params, sw[i]);
            krn[i] = MaterialLaw::krn(params, sw[i]);
            krMin = std::min(krMin, std::min(krw[i], krn[i]));
            krMax = std::max(krMax, std::max(krw[i], krn[i]));
        }

        unsigned int windowNumber = 5;
        gnuplot_.setXRange(lowerSat, upperSat, windowNumber);
        gnuplot_.setYRange(krMin, krMax, windowNumber);
        gnuplot_.setXlabel("wetting phase saturation [-]", windowNumber);
        gnuplot_.setYlabel("relative permeability [-]", windowNumber);
        gnuplot_.addDataSetToPlot(sw, krw, plotName + "_krw", windowNumber);
        gnuplot_.addDataSetToPlot(sw, krn, plotName + "_krn", windowNumber);
        gnuplot_.plot("kr", windowNumber, interaction);
    }

    /*!
     * \brief Plot the gradient of the relative permeabilities
     *
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotdkrdsw(const MaterialLawParams &params,
                    Scalar lowerSat = 0.0,
                    Scalar upperSat = 1.0,
                    std::string plotName = "",
                    bool interaction = true)
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> dkrw_dsw(numIntervals_+1);
        std::vector<Scalar> dkrn_dsw(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar dkrdswMin = 1e100;
        Scalar dkrdswMax = -1e100;
        checkEffectiveSaturation(params, lowerSat, upperSat, plotName + "_dkr_dsw");

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            dkrw_dsw[i] = MaterialLaw::dkrw_dsw(params, sw[i]);
            dkrn_dsw[i] = MaterialLaw::dkrn_dsw(params, sw[i]);
            dkrdswMin = std::min(dkrdswMin, std::min(dkrw_dsw[i], dkrn_dsw[i]));
            dkrdswMax = std::max(dkrdswMax, std::max(dkrw_dsw[i], dkrn_dsw[i]));
        }

        unsigned int windowNumber = 6;
        gnuplot_.setXRange(lowerSat, upperSat, windowNumber);
        gnuplot_.setYRange(dkrdswMin, dkrdswMax, windowNumber);
        gnuplot_.setXlabel("wetting phase saturation [-]", windowNumber);
        gnuplot_.setYlabel("gradient of the kr-Sw function [-]", windowNumber);
        gnuplot_.addDataSetToPlot(sw, dkrw_dsw, plotName + "_dkrw_dsw", windowNumber);
        gnuplot_.addDataSetToPlot(sw, dkrn_dsw, plotName + "_dkrn_dsw", windowNumber);
        gnuplot_.plot("dkrndsw", windowNumber, interaction);
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
    GnuplotInterface<Scalar> gnuplot_;
    int numIntervals_;
};
} // end of namespace

#endif // DUMUX_PLOT_FLUID_MATRIX_LAW_HH
