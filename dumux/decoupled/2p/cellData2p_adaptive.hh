/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_ELEMENTDATA2P_ADAPTIVE_HH
#define DUMUX_ELEMENTDATA2P_ADAPTIVE_HH

#include "cellData2p.hh"

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations for one element
 * @author Markus Wolff
 */

namespace Dumux
{
/*!
 * \ingroup IMPES
 */
//! Class including the variables and data of discretized data of the constitutive relations for one element.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * @tparam TypeTag The Type Tag
 1*/
template<class TypeTag, bool enableCompressibility = GET_PROP_VALUE(TypeTag, EnableCompressibility)>
class CellData2PAdaptive: public CellData2P<TypeTag, enableCompressibility>
{
private:
    typedef CellData2P<TypeTag, enableCompressibility> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef FluxData2P<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

private:
    bool isFront_;

public:

    struct AdaptedValues
    {
        Scalar saturationW;
        Scalar saturationNW;
        Scalar pressW;
        Scalar pressNW;
        Scalar volCorr;
        int count;
        bool front;
        AdaptedValues()
        {
            saturationW = 0.;
            saturationNW = 0.;
            pressW = 0.;
            pressNW = 0.;
            count = 0;
            volCorr = 0;
            front = false;
        }
    };

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    CellData2PAdaptive():
        isFront_(false)
    {}

    bool isFront() const
    {
        return isFront_;
    }

    bool isFront()
    {
        return isFront_;
    }

    void setFront()
    {
        isFront_ = true;
    }

    void getAdaptionValues(AdaptedValues& adaptedValues, const Problem& problem)
    {
        adaptedValues.saturationW = this->saturation(wPhaseIdx);
        adaptedValues.saturationNW = this->saturation(nPhaseIdx);
        adaptedValues.pressW = this->pressure(wPhaseIdx);
        adaptedValues.pressNW = this->pressure(nPhaseIdx);
        adaptedValues.volCorr = this->volumeCorrection();
        adaptedValues.front = isFront_;
    }

    static void getAdaptionValues(AdaptedValues& adaptedValues, AdaptedValues& adaptedValuesFather, const Problem& problem)
    {
        adaptedValuesFather.saturationW += adaptedValues.saturationW / adaptedValues.count;
        adaptedValuesFather.saturationNW += adaptedValues.saturationNW / adaptedValues.count;
        adaptedValuesFather.pressW += adaptedValues.pressW / adaptedValues.count;
        adaptedValuesFather.pressNW += adaptedValues.pressNW / adaptedValues.count;
        adaptedValuesFather.volCorr += adaptedValues.volCorr / adaptedValues.count;
    }

    void setAdaptionValues(AdaptedValues& adaptedValues, const Problem& problem)
    {
        this->setSaturation(wPhaseIdx, adaptedValues.saturationW / adaptedValues.count);
        this->setSaturation(nPhaseIdx, adaptedValues.saturationNW / adaptedValues.count);
        this->setPressure(wPhaseIdx, adaptedValues.pressW / adaptedValues.count);
        this->setPressure(nPhaseIdx, adaptedValues.pressNW / adaptedValues.count);
        this->setUpdate(adaptedValues.volCorr / adaptedValues.count);
        isFront_ = adaptedValues.front;
    }

    static void setAdaptionValues(AdaptedValues& adaptedValues, AdaptedValues& adaptedValuesFather, const Problem& problem)
    {
        adaptedValues.saturationW = adaptedValuesFather.saturationW / adaptedValuesFather.count;
        adaptedValues.saturationNW = adaptedValuesFather.saturationNW / adaptedValuesFather.count;
        adaptedValues.pressW = adaptedValuesFather.pressW / adaptedValuesFather.count;
        adaptedValues.pressNW = adaptedValuesFather.pressNW / adaptedValuesFather.count;
        adaptedValues.volCorr = adaptedValuesFather.volCorr / adaptedValuesFather.count;
    }

};

}
#endif
