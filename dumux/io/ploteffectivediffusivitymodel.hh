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
    PlotEffectiveDiffusivityModel()
    : numIntervals_(1000)
    { }

    /*!
     * \brief Plot the effective diffusion factor-saturation curve
     *
     * \param porosity The porosity of the porous medium
     * \param lowerSat Minimum x-value
     * \param upperSat Maximum x-value
     * \param plotName Name of the plotted curve
     * \param interaction Specifies whether a live output via a gnuplot window is wanted
     */
    void plotdeff(Scalar porosity,
                  Scalar lowerSat = 0.0,
                  Scalar upperSat = 1.0,
                  std::string plotName = "D_eff",
                  bool interaction = true)
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

        unsigned int windowNumber = 7;
        gnuplot_.setXRange(lowerSat, upperSat, windowNumber);
        gnuplot_.setYRange(deffMin, deffMax, windowNumber);
        gnuplot_.setXlabel("phase saturation [-]", windowNumber);
        gnuplot_.setYlabel("effective diffusion/molecular diffusion [-]", windowNumber);
        gnuplot_.addDataSetToPlot(sw, deff, plotName, windowNumber);
        gnuplot_.plot("deff", windowNumber, interaction);
    }

private:
    GnuplotInterface<Scalar> gnuplot_;
    int numIntervals_;
};
} // end of namespace

#endif // DUMUX_PLOT_EFFECTIVE_DIFFUSIVITY_MODEL_HH
