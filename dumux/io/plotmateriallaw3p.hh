// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Interface for plotting the three-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_FLUID_MATRIX_LAW_HH
#define DUMUX_PLOT_FLUID_MATRIX_LAW_HH

#include <cmath>
#include <string>
#include <vector>

namespace Dumux {

// forward declaration
template<class Scalar> class GnuplotInterface;

/*!
 * \ingroup InputOutput
 * \brief Interface for plotting the three-phase fluid-matrix-interaction laws
 * TODO: add theta head pressure plot (porosity and density is needed)
 */
template<class Scalar, class MaterialLaw>
class PlotMaterialLaw
{
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    //! Constructor
    PlotMaterialLaw(bool interaction = true)
    : numIntervals_(1000)
    {}

    /*!
     * \brief Plot the capillary pressure-saturation curve for all  phases
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void addpc(GnuplotInterface<Scalar> &gnuplot,
                const MaterialLawParams &params,
                Scalar lowerSat = 0.0,
                Scalar upperSat = 1.0,
                std::string curveTitle = "")
    {
        addpcgw(gnuplot, params, lowerSat, upperSat, curveTitle);
        addpcnw(gnuplot, params, lowerSat, upperSat, curveTitle);
        addpcgn(gnuplot, params, lowerSat, upperSat, curveTitle);
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the water-gas interphase
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void addpcgw(GnuplotInterface<Scalar> &gnuplot,
                  const MaterialLawParams &params,
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
            pcTemp = MaterialLaw::pcgw(params, swTemp);
            using std::max;
            using std::min;
            if (checkValues_(swTemp, pcTemp))
            {
                sw.push_back(swTemp);
                pc.push_back(pcTemp);
                pcMin = min(pcMin, pcTemp);
                pcMax = max(pcMax, pcTemp);
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("capillary pressure [Pa]");
        gnuplot.addDataSetToPlot(sw, pc, curveTitle + "_pcgw-Sw");
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the water-NAPL interface
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void addpcnw(GnuplotInterface<Scalar> &gnuplot,
                  const MaterialLawParams &params,
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
            pcTemp = MaterialLaw::pcnw(params, swTemp);
            using std::max;
            using std::min;
            if (checkValues_(swTemp, pcTemp))
            {
                sw.push_back(swTemp);
                pc.push_back(pcTemp);
                pcMin = min(pcMin, pcTemp);
                pcMax = max(pcMax, pcTemp);
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("capillary pressure [Pa]");
        gnuplot.addDataSetToPlot(sw, pc, curveTitle + "_pcnw-Sw");
    }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the gas-NAPL interface
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void addpcgn(GnuplotInterface<Scalar> &gnuplot,
                  const MaterialLawParams &params,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string curveTitle = "")
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
            using std::max;
            using std::min;
            if (checkValues_(stTemp, pcTemp))
            {
                st.push_back(stTemp);
                pc.push_back(pcTemp);
                pcMin = min(pcMin, pcTemp);
                pcMax = max(pcMax, pcTemp);
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("capillary pressure [Pa]");
        gnuplot.addDataSetToPlot(st, pc, curveTitle + "_pcgn-St");
    }


    /*!
     * \brief Plot the relative permeabilities
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void addkr(GnuplotInterface<Scalar> &gnuplot,
                const MaterialLawParams &params,
                Scalar lowerSat = 0.0,
                Scalar upperSat = 1.0,
                std::string curveTitle = "")
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
            krwTemp = MaterialLaw::krw(params, swTemp, 0.0);
            krnTemp = MaterialLaw::krn(params, swTemp, 1.0 - swTemp);
            krgTemp = MaterialLaw::krg(params, swTemp, 0.0);
            using std::max;
            using std::min;
            if (checkValues_(swTemp, krwTemp)
                && checkValues_(swTemp, krnTemp)
                && checkValues_(swTemp, krgTemp))
            {
                sw.push_back(swTemp);
                krw.push_back(krwTemp);
                krn.push_back(krnTemp);
                krg.push_back(krgTemp);
                krMin = min({krMin, krwTemp, krnTemp, krgTemp});
                krMax = max({krMax, krwTemp, krnTemp, krgTemp});
            }
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("relative permeability [-]");
        gnuplot.addDataSetToPlot(sw, krw, curveTitle + "_krw");
        gnuplot.addDataSetToPlot(sw, krn, curveTitle + "_krn");
        gnuplot.addDataSetToPlot(sw, krg, curveTitle + "_krg");
    }

    /*!
     * \brief Plot the transition (2P/3P) function
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param curveTitle Name of the plotted curve
     */
    void addPcAlpha(GnuplotInterface<Scalar> &gnuplot,
                     const MaterialLawParams &params,
                     Scalar lowerSat = 0.0,
                     Scalar upperSat = 1.0,
                     std::string curveTitle = "")
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
            using std::max;
            using std::min;
            if (checkValues_(snTemp, alphaTemp))
            {
                sn.push_back(snTemp);
                alpha.push_back(alphaTemp);
                alphaMin = min(alphaMin, alphaTemp);
                alphaMax = max(alphaMax, alphaTemp);
            }
        }

        gnuplot.setXlabel("nonwetting phase saturation [-]");
        gnuplot.setYlabel("transition function [-]");
        gnuplot.addDataSetToPlot(sn, alpha, curveTitle + "_alpha");
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
} // end of namespace Dumux

#endif
