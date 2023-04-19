// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Interface for plotting the multi-component-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_EFFECTIVE_DIFFUSIVITY_MODEL_HH
#define DUMUX_PLOT_EFFECTIVE_DIFFUSIVITY_MODEL_HH

#include <string>
#include <vector>

namespace Dumux {

// forward declaration
template<class Scalar> class GnuplotInterface;

/*!
 * \ingroup InputOutput
 * \brief Interface for plotting the multi-component-matrix-interaction laws
 */
template<class Scalar, class EffectiveDiffusivityModel>
class PlotEffectiveDiffusivityModel
{
public:
    //! Constructor
    PlotEffectiveDiffusivityModel()
    : numIntervals_(1000)
    { }

    /*!
     * \brief Add a effective diffusion factor-saturation data set to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param porosity The porosity
     * \param diffCoeff The binary diffusion coefficient
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void adddeffcurve(GnuplotInterface<Scalar> &gnuplot,
                      Scalar porosity,
                      Scalar diffCoeff,
                      Scalar lowerSat = 0.0,
                      Scalar upperSat = 1.0,
                      std::string curveName = "deff",
                      std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> deff(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            VolumeVariables volVars(sw[i], porosity, diffCoeff);
            deff[i] = EffectiveDiffusivityModel::effectiveDiffusionCoefficient(volVars, 0, 0, 1);
        }

        gnuplot.addDataSetToPlot(sw, deff, curveName, curveOptions);
    }

    //! for point check
    Scalar getEffectiveDiffusionCoefficient(Scalar saturation,
                                            Scalar porosity,
                                            Scalar diffCoeff) const
    {
        VolumeVariables volVars(saturation, porosity, diffCoeff);
        return EffectiveDiffusivityModel::effectiveDiffusionCoefficient(volVars, 0, 0, 1);
    }

private:

    class VolumeVariables
    {
    public:
        VolumeVariables(Scalar saturation, Scalar porosity, Scalar diffCoeff)
        : saturation_(saturation)
        , porosity_(porosity)
        , diffCoeff_(diffCoeff)
        {}

        Scalar saturation(int phaseIdx = 0) const
        { return saturation_; }

        Scalar porosity() const
        { return porosity_; }

        Scalar diffusionCoefficient(const int phaseIdx, const int compIdxI, const int compIdxJ) const
        { return diffCoeff_;}

    private:
        Scalar saturation_;
        Scalar porosity_;
        Scalar diffCoeff_;
    };

    int numIntervals_;
};

} // end namespace Dumux

#endif
