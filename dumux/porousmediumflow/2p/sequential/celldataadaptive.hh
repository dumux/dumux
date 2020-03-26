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
 * \ingroup SequentialTwoPModel
 * \brief  Class including the data of a grid cell needed if an adaptive grid is used.
 */
#ifndef DUMUX_ELEMENTDATA2P_ADAPTIVE_HH
#define DUMUX_ELEMENTDATA2P_ADAPTIVE_HH

#include <dune/grid/utility/persistentcontainer.hh>
#include "celldata.hh"

namespace Dumux {
/*!
 * \ingroup SequentialTwoPModel
 * \brief Class including the data of a grid cell needed if an adaptive grid is used.
 *
 * The class provides model-specific functions needed to adapt the stored cell data to a new (adapted) grid.
 * Additionally, it provides the storage-infrastructure for explicit front tracking.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization for case of compressible flow (<tt>true</tt>) or incompressible flow (<tt>false</tt>)
 */
template<class TypeTag, bool enableCompressibility = getPropValue<TypeTag, Properties::EnableCompressibility>()>
class CellData2PAdaptive: public CellData2P<TypeTag, enableCompressibility>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

public:
    /*!
     * \brief Collection of variables that have to be mapped if the grid is adapted
     * For an immiscible two-phase model, the following data has to be
     * transferred to a new grid.
     */
    struct AdaptedValues
    {
        Scalar saturationW;
        Scalar saturationNw;
        Scalar pressW;
        Scalar pressNw;
        Scalar potW;
        Scalar potNw;
        Scalar volCorr;
        int count;
        AdaptedValues()
        {
            saturationW = 0.;
            saturationNw = 0.;
            pressW = 0.;
            pressNw = 0.;
            potW = 0.;
            potNw = 0.;
            count = 0;
            volCorr = 0;
        }
    };

    using LoadBalanceData = AdaptedValues;

    //! Constructs an adaptive CellData object
    CellData2PAdaptive()
    {}

    /*!
     * \brief Stores values to be adapted in an adaptedValues container
     *
     * Stores values to be adapted from the current CellData objects into
     * the adaptation container in order to be mapped on a new grid.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element to be stored
     */
    void storeAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        adaptedValues.saturationW = this->saturation(wPhaseIdx);
        adaptedValues.saturationNw = this->saturation(nPhaseIdx);
        adaptedValues.pressW = this->pressure(wPhaseIdx);
        adaptedValues.pressNw = this->pressure(nPhaseIdx);
        adaptedValues.potW = this->potential(wPhaseIdx);
        adaptedValues.potNw = this->potential(nPhaseIdx);
        adaptedValues.volCorr = this->volumeCorrection();
    }

    /*!
     * \brief Stores sons entries into father element for averaging
     *
     * Sum up the adaptedValues (sons values) into father element. We store from leaf
     * upwards, so sons are stored first, then cells on the next leaf (=fathers)
     * can be averaged.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param adaptedValuesFather Values to be adapted of father cell
     * \param fatherElement The element of the father
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather,
                                    const Element& fatherElement)
    {
        adaptedValuesFather.saturationW += adaptedValues.saturationW / adaptedValues.count;
        adaptedValuesFather.saturationNw += adaptedValues.saturationNw / adaptedValues.count;
        adaptedValuesFather.pressW += adaptedValues.pressW / adaptedValues.count;
        adaptedValuesFather.pressNw += adaptedValues.pressNw / adaptedValues.count;
        adaptedValuesFather.potW += adaptedValues.potW / adaptedValues.count;
        adaptedValuesFather.potNw += adaptedValues.potNw / adaptedValues.count;
        adaptedValuesFather.volCorr += adaptedValues.volCorr / adaptedValues.count;
    }

    /*!
     * \brief Set adapted values in CellData
     *
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * sequential models.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element where things are stored.
     */
    void setAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        this->setSaturation(wPhaseIdx, adaptedValues.saturationW / adaptedValues.count);
        this->setSaturation(nPhaseIdx, adaptedValues.saturationNw / adaptedValues.count);
        this->setPressure(wPhaseIdx, adaptedValues.pressW / adaptedValues.count);
        this->setPressure(nPhaseIdx, adaptedValues.pressNw / adaptedValues.count);
        this->setPotential(wPhaseIdx, adaptedValues.potW / adaptedValues.count);
        this->setPotential(nPhaseIdx, adaptedValues.potNw / adaptedValues.count);
        this->setUpdate(adaptedValues.volCorr / adaptedValues.count);
    }

    /*!
     * \brief Reconstructs sons entries from data of father cell
     *
     * Reconstructs a new solution from a father cell into a newly
     * generated son cell. New cell is stored into the global
     * adaptationMap.
     *
     * \param adaptionMap Global map storing all values to be adapted
     * \param father Entity Pointer to the father cell
     * \param son Entity Pointer to the newly created son cell
     * \param problem The problem
     */
    static void reconstructAdaptionValues(Dune::PersistentContainer<Grid, AdaptedValues>& adaptionMap,
            const Element& father, const Element& son, const Problem& problem)
    {
        AdaptedValues& adaptedValues = adaptionMap[son];
        AdaptedValues& adaptedValuesFather = adaptionMap[father];
        adaptedValues.saturationW = adaptedValuesFather.saturationW / adaptedValuesFather.count;
        adaptedValues.saturationNw = adaptedValuesFather.saturationNw / adaptedValuesFather.count;
        adaptedValues.pressW = adaptedValuesFather.pressW / adaptedValuesFather.count;
        adaptedValues.pressNw = adaptedValuesFather.pressNw / adaptedValuesFather.count;
        adaptedValues.potW = adaptedValuesFather.potW / adaptedValuesFather.count;
        adaptedValues.potNw = adaptedValuesFather.potNw / adaptedValuesFather.count;
        adaptedValues.volCorr = adaptedValuesFather.volCorr / adaptedValuesFather.count;
    }

};

} // end namespace Dumux
#endif
