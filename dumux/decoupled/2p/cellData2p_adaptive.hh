// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef FluxData2P<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef Dune::FieldVector<Scalar, GridView::dimensionworld> GlobalPosition;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

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
    //! Collection of variables that have to be mapped if the grid is adapted
    /**
     * For an non-compositional two-phase model, the following data has to be
     * transferred to a new grid.
     */
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

    //! Constructs an adaptive CellData object
    /**
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

    //! Stores values to be adapted in an adaptedValues container
    /**
     * Stores values to be adapted from the current CellData objects into
     * the adaptation container in order to be mapped on a new grid.
     *
     * @param adaptedValues Container for model-specific values to be adapted
     * @param problem The problem
     */
    void storeAdaptionValues(AdaptedValues& adaptedValues, const Problem& problem)
    {
        adaptedValues.saturationW = this->saturation(wPhaseIdx);
        adaptedValues.saturationNW = this->saturation(nPhaseIdx);
        adaptedValues.pressW = this->pressure(wPhaseIdx);
        adaptedValues.pressNW = this->pressure(nPhaseIdx);
        adaptedValues.volCorr = this->volumeCorrection();
        adaptedValues.front = isFront_;
    }
    //! Stores sons entries into father element for averageing
    /**
     * Sum up the adaptedValues (sons values) into father element. We store from leaf
     * upwards, so sons are stored first, then cells on the next leaf (=fathers)
     * can be averaged.
     *
     * @param adaptedValues Container for model-specific values to be adapted
     * @param adaptedValuesFather Values to be adapted of father cell
     * @param problem The problem
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues, AdaptedValues& adaptedValuesFather, const Problem& problem)
    {
        adaptedValuesFather.saturationW += adaptedValues.saturationW / adaptedValues.count;
        adaptedValuesFather.saturationNW += adaptedValues.saturationNW / adaptedValues.count;
        adaptedValuesFather.pressW += adaptedValues.pressW / adaptedValues.count;
        adaptedValuesFather.pressNW += adaptedValues.pressNW / adaptedValues.count;
        adaptedValuesFather.volCorr += adaptedValues.volCorr / adaptedValues.count;
    }
    //! Set adapted values in CellData
    /**
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * decoupled models.
     *
     * @param adaptedValues Container for model-specific values to be adapted
     * @param problem The problem
     */
    void setAdaptionValues(AdaptedValues& adaptedValues, const Problem& problem)
    {
        this->setSaturation(wPhaseIdx, adaptedValues.saturationW / adaptedValues.count);
        this->setSaturation(nPhaseIdx, adaptedValues.saturationNW / adaptedValues.count);
        this->setPressure(wPhaseIdx, adaptedValues.pressW / adaptedValues.count);
        this->setPressure(nPhaseIdx, adaptedValues.pressNW / adaptedValues.count);
        this->setUpdate(adaptedValues.volCorr / adaptedValues.count);
        isFront_ = adaptedValues.front;
    }

    //! Reconstructs sons entries from data of father cell
    /**
     * Reconstructs an new solution from a father cell into for a newly
     * generated son cell. New cell is stored into the global
     * adaptationMap.
     *
     * @param adaptionMap Global map storing all values to be adapted
     * @param father Entity Pointer to the father cell
     * @param son Entity Pointer to the newly created son cell
     * @param problem The problem
     */
    static void reconstructAdaptionValues(Dune::PersistentContainer<Grid, AdaptedValues>& adaptionMap,
            const Element& father, const Element& son, const Problem& problem)
    {
        AdaptedValues& adaptedValues = adaptionMap[son];
        AdaptedValues& adaptedValuesFather = adaptionMap[father];
        adaptedValues.saturationW = adaptedValuesFather.saturationW / adaptedValuesFather.count;
        adaptedValues.saturationNW = adaptedValuesFather.saturationNW / adaptedValuesFather.count;
        adaptedValues.pressW = adaptedValuesFather.pressW / adaptedValuesFather.count;
        adaptedValues.pressNW = adaptedValuesFather.pressNW / adaptedValuesFather.count;
        adaptedValues.volCorr = adaptedValuesFather.volCorr / adaptedValuesFather.count;
    }

};

}
#endif
