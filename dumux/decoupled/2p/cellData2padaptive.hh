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
 * @brief  Class including the data of a grid cell needed if an adaptive grid is used.
 */

namespace Dumux
{

/*!
 * \ingroup IMPES
 */
//! Class including the data of a grid cell needed if an adaptive grid is used.
/*! The class provides model-specific functions needed to adapt the stored cell data to a new (adapted) grid.
 * Additionally, it provides the storage-infrastructure for explicit front tracking.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization for case of compressible flow (<tt>true</tt>) or incompressible flow (<tt>false</tt>)
 */
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
    Scalar dt_;

public:
    //! Collection of variables that have to be mapped if the grid is adapted
    /**
     * For an immiscible two-phase model, the following data has to be
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
        AdaptedValues()
        {
            saturationW = 0.;
            saturationNW = 0.;
            pressW = 0.;
            pressNW = 0.;
            count = 0;
            volCorr = 0;
        }
    };

    //! Constructs an adaptive CellData object
    CellData2PAdaptive():
        isFront_(false), dt_(0.0)
    {}

    /*! \brief Track the front
     *
     * Returns true if the cell is located at a fluid-fluid displacement-front
     */
    bool isFront() const
    {
        return isFront_;
    }

    /*! \brief Track the front
     *
     * Returns true if the cell is located at a fluid-fluid displacement-front
     */
    bool isFront()
    {
        return isFront_;
    }

    /*! \brief Reset the front marker
     *
     * Sets front marker to <tt>false</tt>;
     */
    void resetFrontMarker()
    {
        isFront_ = false;
    }

    //! Marks the cell as fluid-fluid displacement-front cell
    void setFront()
    {
        isFront_ = true;
    }

    Scalar dt()
    {
        return dt_;
    }

    Scalar dt() const
    {
        return dt_;
    }

    void setDt(Scalar dt)
    {
        dt_ = dt;
    }

    //! Stores values to be adapted in an adaptedValues container
    /**
     * Stores values to be adapted from the current CellData objects into
     * the adaptation container in order to be mapped on a new grid.
     *
     * @param adaptedValues Container for model-specific values to be adapted
     * @param element The element to be stored
     */
    void storeAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        adaptedValues.saturationW = this->saturation(wPhaseIdx);
        adaptedValues.saturationNW = this->saturation(nPhaseIdx);
        adaptedValues.pressW = this->pressure(wPhaseIdx);
        adaptedValues.pressNW = this->pressure(nPhaseIdx);
        adaptedValues.volCorr = this->volumeCorrection();
    }
    //! Stores sons entries into father element for averageing
    /**
     * Sum up the adaptedValues (sons values) into father element. We store from leaf
     * upwards, so sons are stored first, then cells on the next leaf (=fathers)
     * can be averaged.
     *
     * @param adaptedValues Container for model-specific values to be adapted
     * @param adaptedValuesFather Values to be adapted of father cell
     * @param fatherElement The element of the father
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather,
                                    const Element& fatherElement)
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
     * @param element The element where things are stored.
     */
    void setAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        this->setSaturation(wPhaseIdx, adaptedValues.saturationW / adaptedValues.count);
        this->setSaturation(nPhaseIdx, adaptedValues.saturationNW / adaptedValues.count);
        this->setPressure(wPhaseIdx, adaptedValues.pressW / adaptedValues.count);
        this->setPressure(nPhaseIdx, adaptedValues.pressNW / adaptedValues.count);
        this->setUpdate(adaptedValues.volCorr / adaptedValues.count);
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
