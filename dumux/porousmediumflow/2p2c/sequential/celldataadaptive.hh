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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Class including the variables and data of discretized data of the constitutive relations for one element
 */
#ifndef DUMUX_ELEMENTDATA2P2C_ADAPTIVE_HH
#define DUMUX_ELEMENTDATA2P2C_ADAPTIVE_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/porousmediumflow/2p2c/sequential/celldata.hh>
#include <dumux/porousmediumflow/2p2c/sequential/celldatamultiphysics.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Class including the data of a grid cell needed if an adaptive grid is used.
 *
 * The class provides model-specific functions needed to adapt the stored cell data to a new (adapted) grid.
 * Additionally, it provides the storage-infrastructure for explicit front tracking.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CellData2P2CAdaptive: public CellData2P2CMultiPhysics<TypeTag>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    enum
    {
        dim = GridView::dimension
    };

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx, nCompIdx = Indices::nCompIdx
    };
    enum
    {
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    //! gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static constexpr int pressureType = getPropValue<TypeTag, Properties::PressureFormulation>();
    int upwindError_[numPhases];

public:
    /*!
     * \brief A container for all necessary variables to map an old solution to a new grid
     * If the primary variables (pressure, total concentrations) are mapped to a new grid,
     * the secondary variables can be calulated. For the mapping between sons and father, it
     * is in addition necessary to know about how many sons live in each father ("count").
     * While only one phase pressure is PV, the method updateMaterialLaws() that
     * recalculates the secondary variables needs both phase pressures (initiated via the
     * capillary pressure of the last time step) to start the iteration to determine
     * appropriate phase composition and pc. Hence, both phase pressures are mapped
     * to the new solution.
     */
    struct AdaptedValues
    {
        Dune::FieldVector<Scalar,2> totalConcentration_; //!< Transport primary variables
        Dune::FieldVector<Scalar,2> pressure_; //!< Pressure primary variables
        //! Storage for volume derivatives, as transport estimate on old pressure field is not correct after refinement
        Dune::FieldVector<Scalar,3> volumeDerivatives_;
        Scalar cellVolume; //!< Cell volume for transformation of volume-specific primary variables
        FluxData2P2C<TypeTag> fluxData_; //!< Information on old flux direction
        int subdomain; //!< subdomain
        int count; //!< counts the number of cells averaged
        int isRefined; //!< Indicates if cell is refined
        AdaptedValues()
        {
            totalConcentration_=0.;
            pressure_ = 0.;
            volumeDerivatives_ =0.;
            cellVolume = 0.;
            subdomain=-1;
            count = 0;
            isRefined = false;
        }
    };


    //! Constructs an adaptive CellData object
    CellData2P2CAdaptive() : CellData2P2CMultiPhysics<TypeTag>()
    {
        for (int i = 0; i < numPhases;i++)
            upwindError_[i] = 0;
    }

    /*!
     * \brief Stores leaf cell primary variables to transfer to new indexing
     * Stores values to be adapted from the current CellData objects into
     * the adaptation container in order to be mapped on a new grid.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element to be stored
     */
    void storeAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        adaptedValues.totalConcentration_[wCompIdx]= this->massConcentration(wCompIdx);
        adaptedValues.totalConcentration_[nCompIdx]= this->massConcentration(nCompIdx);
        adaptedValues.pressure_[wPhaseIdx] = this->pressure(wPhaseIdx);
        adaptedValues.pressure_[nPhaseIdx] = this->pressure(nPhaseIdx);
        adaptedValues.volumeDerivatives_[wCompIdx] = this->dv(wPhaseIdx);
        adaptedValues.volumeDerivatives_[nCompIdx] = this->dv(nPhaseIdx);
        adaptedValues.volumeDerivatives_[2] = this->dv_dp();
        adaptedValues.cellVolume = element.geometry().volume();
        adaptedValues.subdomain = this->subdomain();
        adaptedValues.fluxData_=this->fluxData();
    }

    /*!
     * \brief Adds cell information to father element for possible averaging / coarsening
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
        adaptedValuesFather.totalConcentration_[wCompIdx]
                += adaptedValues.cellVolume* adaptedValues.totalConcentration_[wCompIdx];
        adaptedValuesFather.totalConcentration_[nCompIdx]
                += adaptedValues.cellVolume* adaptedValues.totalConcentration_[nCompIdx];
        // if all cells are summed up, re-convert mass into total concentrations
        Scalar fatherVolume = fatherElement.geometry().volume();
        if(adaptedValuesFather.count == 1 << dim)
        {
            adaptedValuesFather.totalConcentration_[wCompIdx] /= fatherVolume;
            adaptedValuesFather.totalConcentration_[nCompIdx] /= fatherVolume;
        }
        adaptedValuesFather.cellVolume = fatherVolume;

        adaptedValuesFather.pressure_[wPhaseIdx] += adaptedValues.pressure_[wPhaseIdx] / adaptedValues.count;
        adaptedValuesFather.pressure_[nPhaseIdx] += adaptedValues.pressure_[nPhaseIdx] / adaptedValues.count;
        // apply maximum complexity for new cell
        using std::max;
        adaptedValuesFather.subdomain = max(adaptedValuesFather.subdomain, adaptedValues.subdomain);
    }

    /*!
     * \brief Set adapted values in CellData
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * sequential models.
     * In new cells, update estimate does not give meaningful results. We therefore
     * copy volume derivatives from old time step, and indicate that those are already availabe.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element where things are stored.
     */
    void setAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        // in new cells, there is a cellData object, but not yet a fluidstate.
        // To write the adapted values in this fluidstate, we have to enshure that
        // it was created. We therefore specify the type of FluidState that should be stored.
        this->setSubdomainAndFluidStateType(adaptedValues.subdomain);
        this->setMassConcentration(wCompIdx,
                adaptedValues.totalConcentration_[wCompIdx]);
        this->setMassConcentration(nCompIdx,
                adaptedValues.totalConcentration_[nCompIdx]);
        this->setPressure(wPhaseIdx,
                adaptedValues.pressure_[wPhaseIdx] / adaptedValues.count);
        this->setPressure(nPhaseIdx,
                adaptedValues.pressure_[nPhaseIdx] / adaptedValues.count);

        //copy flux directions
        this->fluxData()=adaptedValues.fluxData_;

        //indicate if cell was just refined in this TS
        this->wasRefined()= adaptedValues.isRefined;
        if(adaptedValues.isRefined)
        {
            this->dv(wPhaseIdx) = adaptedValues.volumeDerivatives_[wCompIdx];
            this->dv(nPhaseIdx) = adaptedValues.volumeDerivatives_[nCompIdx];
            this->dv_dp() = adaptedValues.volumeDerivatives_[2];
            this->volumeDerivativesAvailable(true);
        }
        else
            this->volumeDerivativesAvailable(false);  // recalculate volume derivatives
    }

    /*!
     * \brief Reconstructs sons entries from data of father cell
     * Reconstructs an new solution from a father cell into for a newly
     * generated son cell. The new cell is stored into the global
     * adaptationMap.
     *
     * \param adaptationMap Global map storing all values to be adapted
     * \param father Entity Pointer to the father cell
     * \param son Entity Pointer to the newly created son cell
     * \param problem The problem
     */
    static void reconstructAdaptionValues(Dune::PersistentContainer<Grid, AdaptedValues>& adaptationMap,
            const Element& father, const Element& son, const Problem& problem)
    {
        // acess father and son
        AdaptedValues& adaptedValues = adaptationMap[son];
        AdaptedValues& adaptedValuesFather = adaptationMap[father];

        adaptedValues.subdomain = adaptedValuesFather.subdomain;
        adaptedValues.totalConcentration_[wCompIdx]
                   = adaptedValuesFather.totalConcentration_[wCompIdx] / adaptedValuesFather.count;
        adaptedValues.totalConcentration_[nCompIdx]
                   = adaptedValuesFather.totalConcentration_[nCompIdx] / adaptedValuesFather.count;

        // Introduce a hydrostatic pressure distribution inside the father cell
        Scalar gTimesHeight = problem.gravity()
                    * (son.geometry().center() - father.geometry().center());
        gTimesHeight *= (adaptedValuesFather.totalConcentration_[wCompIdx]+adaptedValuesFather.totalConcentration_[nCompIdx])/
                         problem.spatialParams().porosity(son);
//        int globalIdxSon = problem.variables().index(son);


        // recalculate new Primary variable and store pc (i.e. other pressure)
        if(pressureType == wPhaseIdx)
        {
            // recompute pressure and pc
            Scalar pressure = adaptedValuesFather.pressure_[wPhaseIdx] / adaptedValuesFather.count;
            Scalar pc = adaptedValuesFather.pressure_[nPhaseIdx] / adaptedValuesFather.count
                        - pressure;

            adaptedValues.pressure_[wPhaseIdx]
                    = pressure + gTimesHeight  ;
            adaptedValues.pressure_[nPhaseIdx]
                    = pc + adaptedValues.pressure_[wPhaseIdx];
        }
        else
        {
            // recompute pressure and pc
            Scalar pressure = adaptedValuesFather.pressure_[nPhaseIdx] / adaptedValuesFather.count;
            Scalar pc = pressure
                        - adaptedValuesFather.pressure_[wPhaseIdx] / adaptedValuesFather.count;

            adaptedValues.pressure_[nPhaseIdx]
                    = pressure + gTimesHeight  ;
            adaptedValues.pressure_[wPhaseIdx]
                    = adaptedValues.pressure_[nPhaseIdx] - pc;
        }
        // copy volume derivatives and compressibility, and indicate that cell was refined
        adaptedValues.volumeDerivatives_[wCompIdx] = adaptedValuesFather.volumeDerivatives_[wCompIdx];
        adaptedValues.volumeDerivatives_[nCompIdx] = adaptedValuesFather.volumeDerivatives_[nCompIdx];
        adaptedValues.volumeDerivatives_[2] = adaptedValuesFather.volumeDerivatives_[2];
        adaptedValues.isRefined = true;
    }
};

} // end namespace Dumux
#endif
