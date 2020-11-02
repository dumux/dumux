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
 * \ingroup InputOutput
 * \brief Interface for plotting the two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_FLUID_MATRIX_LAW_HH
#define DUMUX_PLOT_FLUID_MATRIX_LAW_HH

#warning "This header is deprecated and will be removed after 3.3. Use new 2p material laws and plot tools from io/plotpckrsw.hh"

#include <cmath>
#include <vector>
#include <string>

namespace Dumux {

// forward declaration
template<class Scalar> class GnuplotInterface;

/*!
 * \ingroup InputOutput
 * \brief Interface for plotting the two-phase fluid-matrix-interaction laws
 */
template<class Scalar, class MaterialLaw>
class PlotMaterialLaw
{
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    //! Constructor
    PlotMaterialLaw()
    : numIntervals_(1000)
    { }

    /*!
     * \brief Add a capillary pressure-saturation data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addpcswcurve(GnuplotInterface<Scalar> &gnuplot,
                      const MaterialLawParams &params,
                      Scalar lowerSat = 0.0,
                      Scalar upperSat = 1.0,
                      std::string curveName = "pc-Sw",
                      std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            pc[i] = MaterialLaw::pc(params, sw[i]);
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("capillary pressure [Pa]");
        gnuplot.addDataSetToPlot(sw, pc, curveName, curveOptions);
    }

    /*!
     * \brief Add a capillary pressure-saturation data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addLog10PcSwCurve(GnuplotInterface<Scalar> &gnuplot,
                           const MaterialLawParams &params,
                           Scalar lowerSat = 0.0,
                           Scalar upperSat = 1.0,
                           std::string curveName = "log10_pc-Sw.dat",
                           std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> log10pc(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            const Scalar pc = std::log10(MaterialLaw::pc(params, sw[i]));
            log10pc[i] = std::isnan(pc) ? 0.0 : pc;
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("log10 of capillary pressure [Pa]");
        gnuplot.addDataSetToPlot(sw, log10pc, curveName, curveOptions);
    }

    /*!
     * \brief Add a saturation-capillary pressure data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerpc Minimum x-value for data set
     * \param upperpc Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addswpccurve(GnuplotInterface<Scalar> &gnuplot,
                      const MaterialLawParams &params,
                      Scalar lowerpc = 0.0,
                      Scalar upperpc = 5000.0,
                      std::string curveName = "Sw-pc",
                      std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> pc;
        Scalar pcInterval = upperpc - lowerpc;

        Scalar pcTemp, swTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            pcTemp = lowerpc + pcInterval * Scalar(i) / Scalar(numIntervals_);
            swTemp = MaterialLaw::sw(params, pcTemp);
            if (checkValues_(pcTemp, swTemp))
            {
                pc.push_back(pcTemp);
                sw.push_back(swTemp);
            }
        }

        gnuplot.setXlabel("capillary pressure [Pa]");
        gnuplot.setYlabel("wetting phase saturation [-]");
        gnuplot.addDataSetToPlot(pc, sw, curveName, curveOptions);
    }

    /*!
     * \brief Add a saturation-capillary pressure data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerpc Minimum x-value for data set
     * \param upperpc Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addSwLog10PcCurve(GnuplotInterface<Scalar> &gnuplot,
                           const MaterialLawParams &params,
                           Scalar lowerpc = 1,
                           Scalar upperpc = 1e9,
                           std::string curveName = "Sw-pc",
                           std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> pc;
        Scalar pcInterval = std::log10(upperpc) - std::log10(lowerpc);

        Scalar pcTemp, swTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            pcTemp = std::log10(lowerpc) + pcInterval * Scalar(i) / Scalar(numIntervals_);
            swTemp = MaterialLaw::sw(params, std::pow(10, pcTemp));
            pc.push_back(pcTemp);
            sw.push_back(swTemp);
        }

        gnuplot.setXlabel("log10 capillary pressure [Pa]");
        gnuplot.setYlabel("wetting phase saturation [-]");
        gnuplot.addDataSetToPlot(pc, sw, curveName, curveOptions);
    }

    /*!
     * \brief Add a capillary pressure-saturation gradient data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void adddpcdswcurve(GnuplotInterface<Scalar> &gnuplot,
                        const MaterialLawParams &params,
                        Scalar lowerSat = 0.0,
                        Scalar upperSat = 1.0,
                        std::string curveName = "dpcdsw",
                        std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> dpcdsw;
        Scalar satInterval = upperSat - lowerSat;

        Scalar swTemp, dpcdswTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            swTemp = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            dpcdswTemp = MaterialLaw::dpc_dsw(params, swTemp);
            if (checkValues_(swTemp, dpcdsw))
            {
                sw.push_back(swTemp);
                dpcdsw.push_back(dpcdswTemp);
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("gradient of the pc-Sw curve [Pa]");
        gnuplot.addDataSetToPlot(sw, dpcdsw, curveName, curveOptions);
    }

    /*!
     * \brief Add a saturation-capillary pressure gradient data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerpc Minimum x-value for data set
     * \param upperpc Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void adddswdpccurve(GnuplotInterface<Scalar> &gnuplot,
                        const MaterialLawParams &params,
                        Scalar lowerpc = 0.0,
                        Scalar upperpc = 5000.0,
                        std::string curveName = "dswdpc",
                        std::string curveOptions = "w l")
    {
        std::vector<Scalar> pc;
        std::vector<Scalar> dswdpc;
        Scalar pcInterval = upperpc - lowerpc;

        Scalar dswdpcTemp, pcTemp = 0.0;
        for (int i = 0; i <= numIntervals_; i++)
        {
            pcTemp = lowerpc + pcInterval * Scalar(i) / Scalar(numIntervals_);
            dswdpcTemp = MaterialLaw::dsw_dpc(params, pcTemp);
            if (checkValues_(pcTemp, dswdpcTemp))
            {
                pc.push_back(pcTemp);
                dswdpc.push_back(dswdpcTemp);
            }
        }

        gnuplot.setXlabel("capillary pressure [Pa]");
        gnuplot.setYlabel("gradient of the Sw-pc curve [1/Pa]");
        gnuplot.addDataSetToPlot(pc, dswdpc, curveName, curveOptions);
    }

    /*!
     * \brief Add relative permeabilities data sets to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addkrcurves(GnuplotInterface<Scalar> &gnuplot,
                      const MaterialLawParams &params,
                      Scalar lowerSat = 0.0,
                      Scalar upperSat = 1.0,
                      std::string curveName = "kr",
                      std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> krw;
        std::vector<Scalar> krn;
        Scalar satInterval = upperSat - lowerSat;

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
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("relative permeability [-]");
        gnuplot.addDataSetToPlot(sw, krw, curveName + "_krw", curveOptions);
        gnuplot.addDataSetToPlot(sw, krn, curveName + "_krn", curveOptions);
    }

    /*!
     * \brief Add relative permeabilities gradients data sets to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void adddkrdswcurves(GnuplotInterface<Scalar> &gnuplot,
                         const MaterialLawParams &params,
                         Scalar lowerSat = 0.0,
                         Scalar upperSat = 1.0,
                         std::string curveName = "dkrndsw",
                         std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw;
        std::vector<Scalar> dkrw_dsw;
        std::vector<Scalar> dkrn_dsw;
        Scalar satInterval = upperSat - lowerSat;

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
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("gradient of the kr-Sw function [-]");
        gnuplot.addDataSetToPlot(sw, dkrw_dsw, curveName + "_dkrw_dsw", curveOptions);
        gnuplot.addDataSetToPlot(sw, dkrn_dsw, curveName + "_dkrn_dsw", curveOptions);
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
        using std::isnan;
        using std::isinf;
        return !isnan(value1) && !isinf(value1)
               && !isnan(value2) && !isinf(value2);
    }

    int numIntervals_;
};

} // end namespace Dumux

#endif
