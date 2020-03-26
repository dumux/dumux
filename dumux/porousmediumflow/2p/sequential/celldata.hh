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
 * \brief  Class including data of one grid cell
 */
#ifndef DUMUX_ELEMENTDATA2P_HH
#define DUMUX_ELEMENTDATA2P_HH

#include "properties.hh"
#include "fluxdata.hh"

namespace Dumux {
template<class TypeTag>
class FluxData2P;

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class including data of one grid cell.
 *
 * The variables of two-phase flow, which are phase pressures and saturations are stored in this class.
 * Further, resulting cell values for constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure are stored.
 * Additionally, data assigned to cell-cell interfaces, so-called flux-data are stored.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization for case of compressible flow (<tt>true</tt>) or incompressible flow (<tt>false</tt>)
 */
template<class TypeTag, bool enableCompressibility>
class CellData2P;

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class including the variables and data of discretized data of the constitutive relations for one grid cell.
 *
 * The variables of two-phase flow, which are phase pressures and saturations are stored in this class.
 * Further, resulting cell values for constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure are stored.
 * Additionally, data assigned to cell-cell interfaces, so-called flux-data are stored.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization: in case of incompressible flow bool = <tt>false</tt>
 */
template<class TypeTag>
class CellData2P<TypeTag, false>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluxData = FluxData2P<TypeTag>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };
private:
    Scalar saturation_[numPhases];
    Scalar pressure_[numPhases];
    Scalar potential_[numPhases];
    Scalar mobility_[numPhases];
    Scalar fracFlowFunc_[numPhases];

    Scalar capillaryPressure_;
    Scalar update_;

    FluxData fluxData_;
public:

    //! Constructs a CellData2P object
    CellData2P()
    {
        for (int i = 0; i < numPhases;i++)
        {
            saturation_[i] = 0.0;
            pressure_[i] = 0.0;
            potential_[i] = 0.0;
            mobility_[i] = 0.0;
            fracFlowFunc_[i] = 0.0;
        }
        capillaryPressure_ = 0.0;
        update_ = 0.0;
    }

    //! Returns the flux data of the cell
    FluxData& fluxData()
    {
        return fluxData_;
    }

    //! Returns the flux data of the cell
    const FluxData& fluxData() const
    {
        return fluxData_;
    }

    //! \cond \private
    //! Fluidstates are not stored for the incompressible model, however the function is needed to avoid compiler errors
    FluidState& fluidState()
    {
        DUNE_THROW(Dune::NotImplemented,"fluid states not stored in cell data of incompressible models!");
    }

    //! Fluidstates are not stored for the incompressible model, however the function is needed to avoid compiler errors
    const FluidState& fluidState() const
    {
        DUNE_THROW(Dune::NotImplemented,"fluid states not stored in cell data of incompressible models!");
    }
    //! \endcond

    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar pressure(int phaseIdx)
    {
        return pressure_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar pressure(int phaseIdx) const
    {
        return pressure_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     * \param press Phase pressure which is stored
     */
    void setPressure(int phaseIdx, Scalar press)
    {
        pressure_[phaseIdx] = press;
    }

    //! Returns the global pressure of the cell
    Scalar globalPressure()
    {
        return pressure_[wPhaseIdx];
    }

    //! Returns the global pressure of the cell
    Scalar globalPressure() const
    {
        return pressure_[wPhaseIdx];
    }

    /*!
     * \brief Sets the cell global pressure
     *
     * \param press Global pressure which is stored
     */
    void setGlobalPressure(Scalar press)
    {
        pressure_[wPhaseIdx] = press;
    }

    /*!
     * \brief Returns the cell phase potential
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar potential(int phaseIdx)
    {
        return potential_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase potential
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar potential(int phaseIdx) const
    {
        return potential_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase potential
     *
     * \param phaseIdx Index of a fluid phase
     * \param pot Phase potential which is stored
     */
    void setPotential(int phaseIdx, Scalar pot)
    {
        potential_[phaseIdx] = pot;
    }

    /*!
     * \brief Returns the cell phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar saturation(int phaseIdx)
    {
        return saturation_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar saturation(int phaseIdx) const
    {
        return saturation_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     * \param sat Phase saturation which is stored
     */
    void setSaturation(int phaseIdx, Scalar sat)
    {
        saturation_[phaseIdx] = sat;
    }

    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the cell phase mobility
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar mobility(int phaseIdx)
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase mobility
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase mobility
     *
     * \param phaseIdx Index of a fluid phase
     * \param mobility Phase mobility with which is stored
     */
    void setMobility(int phaseIdx, Scalar mobility)
    {
        mobility_[phaseIdx] = mobility;
    }

    /*!
     * \brief Returns the cell phase fractional flow function
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar fracFlowFunc(int phaseIdx)
    {
        return fracFlowFunc_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase fractional flow function
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar fracFlowFunc(int phaseIdx) const
    {
        return fracFlowFunc_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase fractional flow function
     *
     * \param phaseIdx Index of a fluid phase
     * \param fracFlowFunc Phase fractional flow function which is stored
     */
    void  setFracFlowFunc(int phaseIdx, Scalar fracFlowFunc)
    {
        fracFlowFunc_[phaseIdx] = fracFlowFunc;
    }

    //! Returns the cell capillary pressure
    Scalar capillaryPressure()
    {
        return capillaryPressure_;
    }

    //! Returns the cell capillary pressure
    Scalar capillaryPressure() const
    {
        return capillaryPressure_;
    }

    /*!
     * \brief Sets the cell capillary pressure
     *
     * \param pc Capillary pressure which is stored
     */
    void setCapillaryPressure(Scalar pc)
    {
        capillaryPressure_ = pc;
    }

    /*!
     * \brief Store transport update
     *
     * \param update Transport update of the cell
     */
    void setUpdate(Scalar update)
    {
        update_ = update;
    }

    //! Returns the cell volume correction needed in the pressure equation
    Scalar volumeCorrection()
    {
        return update_;
    }
    //! Returns the cell volume correction needed in the pressure equation
    Scalar volumeCorrection() const
    {
        return update_;
    }

    //! \cond \private
    //densities are not stored for the incompressible model, however the function is needed to avoid compiler errors
    Scalar density(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,"density function not implemented in cell data of incompressible models!");
    }
    Scalar density(int phaseIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,"density function not implemented in cell data of incompressible models!");
    }

    //viscosities are not stored for the incompressible model, however the function is needed to avoid compiler errors
    Scalar viscosity(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,"viscosity function not implemented in cell data of incompressible models!");
    }

    Scalar viscosity(int phaseIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,"viscosity function not implemented in cell data of incompressible models!");
    }
    //! \endcond
};

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class including the variables and data of discretized data of the constitutive relations for one grid cell.
 *
 * The variables of two-phase flow, which are phase pressures and saturations are stored in this class.
 * Further, resulting cell values for constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure are stored.
 * Additionally, data assigned to cell-cell interfaces, so-called flux-data are stored.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization: In case of compressible flow bool = <tt>true</tt>
 */
template<class TypeTag>
class CellData2P<TypeTag, true>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluxData = FluxData2P<TypeTag>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };
private:
    Scalar potential_[numPhases];
    Scalar mobility_[numPhases];
    Scalar fracFlowFunc_[numPhases];
    Scalar update_;

    FluxData fluxData_;
    FluidState fluidState_;
public:

    //! Constructs a CellData2P object
    CellData2P()
    {
        for (int i = 0; i < numPhases;i++)
        {
            mobility_[i] = 0.0;
            potential_[i] = 0.0;
            fracFlowFunc_[i] = 0.0;
        }
        update_ = 0.0;
    }

    //! Returns the flux data of the cell
    FluxData& fluxData()
    {
        return fluxData_;
    }

    //! Returns the flux data of the cell
    const FluxData& fluxData() const
    {
        return fluxData_;
    }

    //! Returns the FluidState object for this cell
    FluidState& fluidState()
    {
        return fluidState_;
    }
    //! Returns the FluidState object for this cell
    const FluidState& fluidState() const
    {
        return fluidState_;
    }

    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar pressure(int phaseIdx)
    {
        return fluidState_.pressure(phaseIdx);
    }

    /*!
     * \brief Returns the cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar pressure(int phaseIdx) const
    {
        return fluidState_.pressure(phaseIdx);
    }

    /*!
     * \brief Sets the cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     * \param press Phase pressure which is stored
     */
    void setPressure(int phaseIdx, Scalar press)
    {
        fluidState_.setPressure(phaseIdx, press);
    }

    /*!
     * \brief Returns the cell phase potential
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar potential(int phaseIdx)
    {
        return potential_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase potential
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar potential(int phaseIdx) const
    {
        return potential_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase potential
     *
     * \param phaseIdx Index of a fluid phase
     * \param pot Phase potential which is stored
     */
    void setPotential(int phaseIdx, Scalar pot)
    {
        potential_[phaseIdx] = pot;
    }

    //! \cond \private
    // global pressure is not supported in the compressible case, however the function hast to be defined to avoid compiler errors!
    Scalar globalPressure()
    {
        DUNE_THROW(Dune::NotImplemented,"no global pressure defined for compressible models!");
    }

    // global pressure is not supported in the compressible case, however the function hast to be defined to avoid compiler errors!
    Scalar globalPressure() const
    {
        DUNE_THROW(Dune::NotImplemented,"no global pressure defined for compressible models!");
    }

    // global pressure is not supported in the compressible case, however the function hast to be defined to avoid compiler errors!
    void setGlobalPressure(Scalar press)
    {
        DUNE_THROW(Dune::NotImplemented,"no global pressure defined for compressible models!");
    }
    //! \endcond

    /*!
     * \brief Returns the cell phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar saturation(int phaseIdx)
    {
        return fluidState_.saturation(phaseIdx);
    }

    /*!
     * \brief Returns the cell phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar saturation(int phaseIdx) const
    {
        return fluidState_.saturation(phaseIdx);
    }

    /*!
     * \brief Sets the cell phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     * \param sat Phase saturation which is stored
     */
    void setSaturation(int phaseIdx, Scalar sat)
    {
        fluidState_.setSaturation(phaseIdx, sat);
    }

    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the cell phase mobility
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar mobility(int phaseIdx)
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase mobility
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase mobility
     *
     * \param phaseIdx Index of a fluid phase
     * \param mobility Phase mobility with which is stored
     */
    void setMobility(int phaseIdx, Scalar mobility)
    {
        mobility_[phaseIdx] = mobility;
    }

    /*!
     * \brief Returns the cell phase fractional flow function
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar fracFlowFunc(int phaseIdx)
    {
        return fracFlowFunc_[phaseIdx];
    }

    /*!
     * \brief Returns the cell phase fractional flow function
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar fracFlowFunc(int phaseIdx) const
    {
        return fracFlowFunc_[phaseIdx];
    }

    /*!
     * \brief Sets the cell phase fractional flow function
     *
     * \param phaseIdx Index of a fluid phase
     * \param fracFlowFunc Phase fractional flow function which is stored
     */
    void  setFracFlowFunc(int phaseIdx, Scalar fracFlowFunc)
    {
        fracFlowFunc_[phaseIdx] = fracFlowFunc;
    }

    //! Returns the cell capillary pressure
    Scalar capillaryPressure()
    {
        return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx);
    }

    //! Returns the cell capillary pressure
    Scalar capillaryPressure() const
    {
        return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx);
    }

    //! \cond \private
    //in the compressible model capillary pressure is not stored explicitly, however the function is needed to avoid compiler errors
    void setCapillaryPressure(Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented,"no capillary pressure stored for compressible models!");
    }
    //! \endcond

    /*!
     * \brief Returns the cell phase density
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar density(int phaseIdx)
    {
        return fluidState_.density(phaseIdx);
    }

    /*!
     * \brief Returns the cell phase density
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar density(int phaseIdx) const
    {
        return fluidState_.density(phaseIdx);
    }

    /*!
     * \brief Sets the cell phase density
     *
     * \param phaseIdx Index of a fluid phase
     * \param density Phase density which is stored
     */
    void setDensity(int phaseIdx, Scalar density)
    {
        fluidState_.setDensity(phaseIdx, density);
    }

    /*!
     * \brief Returns the cell phase viscosity
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar viscosity(int phaseIdx)
    {
        return fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the cell phase viscosity
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar viscosity(int phaseIdx) const
    {
        return fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Sets the cell phase viscosity
     *
     * \param phaseIdx Index of a fluid phase
     * \param viscosity Phase viscosity which is stored
     */
    void setViscosity(int phaseIdx, Scalar viscosity)
    {
        fluidState_.setViscosity(phaseIdx, viscosity);
    }

    /*!
     * \brief Store transport update
     *
     * \param update Transport update of the cell
     */
    void setUpdate(Scalar update)
    {
        update_ = update;
    }

    //! Returns the cell volume correction needed in the pressure equation
    Scalar volumeCorrection()
    {
        return update_;
    }

    //! Returns the cell volume correction needed in the pressure equation
    Scalar volumeCorrection() const
    {
        return update_;
    }

};
} // end namespace Dumux
#endif
