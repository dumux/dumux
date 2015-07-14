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
 * \brief Interface for plotting the multi-component-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_EFFECTIVE_DIFFUSIVITY_MODEL_HH
#define DUMUX_PLOT_EFFECTIVE_DIFFUSIVITY_MODEL_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/io/gnuplotinterface.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(EffectiveDiffusivityModel);
NEW_PROP_TAG(Scalar);
}

/*!
 * \brief Interface for plotting the multi-component-matrix-interaction laws
 */
template<class TypeTag>
class PlotEffectiveDiffusivityModel
{
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    //! Constructor
    PlotEffectiveDiffusivityModel(bool interaction = true)
    : numIntervals_(1000)
    {
        gnuplot_.setInteraction(interaction);
    }

    DUNE_DEPRECATED_MSG("plotdeff() has changed signature")
    void plotdeff(Scalar porosity,
                  Scalar lowerSat,
                  Scalar upperSat,
                  std::string plotName,
                  bool interaction)
    { }

    /*!
     * \brief Plot the effective diffusion factor-saturation curve
     *
     * \param porosity The porosity of the porous medium
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     */
    void plotdeff(Scalar porosity,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> deff(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;
        Scalar deffMin = 1e100;
        Scalar deffMax = -1e100;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            deff[i] = EffectiveDiffusivityModel::effectiveDiffusivity(porosity, sw[i],
                                                                      1.0 /*Diffusion Coefficient*/);
            deffMin = std::min(deffMin, deff[i]);
            deffMax = std::max(deffMax, deff[i]);
        }

        gnuplot_.setXRange(lowerSat, upperSat);
        gnuplot_.setYRange(deffMin, deffMax);
        gnuplot_.setXlabel("phase saturation [-]");
        gnuplot_.setYlabel("effective diffusion/molecular diffusion [-]");
        gnuplot_.addDataSetToPlot(sw, deff, plotName + "_d_eff");
        gnuplot_.plot("deff");
    }

private:
    int numIntervals_;
    GnuplotInterface<Scalar> gnuplot_;
};
} // end of namespace

#endif // DUMUX_PLOT_EFFECTIVE_DIFFUSIVITY_MODEL_HH
