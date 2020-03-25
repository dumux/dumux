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
 * \brief Storage container for discretized data of the constitutive relations for one element
 */

#ifndef DUMUX_ELEMENTDATA2P2C_HH
#define DUMUX_ELEMENTDATA2P2C_HH

#include "fluxdata.hh"

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Storage container for discretized data of the constitutive relations for one element
 *
 * This class stores all cell-centered (FV-Scheme) values for sequential compositional two-phase flow
 * models that are used by both pressure and transport model. All fluid data are already stored in the
 * fluidstate, so the CellData contains the fluidstate object for the current element.
 * At the moment, the compositional model does not use fluxVariables that are stored on the interfaces.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CellData2P2C
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxData = FluxData2P2C<TypeTag>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx
    };
    enum
    {
        numPhases = getPropValue<TypeTag, Properties::NumPhases>(),
        numComponents = getPropValue<TypeTag, Properties::NumComponents>()
    };
protected:
    // primary variable (phase pressure has to be stored in fluidstate)
    Scalar massConcentration_[numComponents];

    Scalar mobility_[numPhases];

    //numerical quantities
    Scalar volumeError_;
    Scalar errorCorrection_;
    Scalar dv_dp_;
    Scalar dv_[numComponents];
    bool volumeDerivativesAvailable_;
    bool wasRefined_;

    int globalIdx_;
    Scalar perimeter_;

    std::shared_ptr<FluidState> fluidState_;
    FluxData fluxData_;
public:

    //! Constructor for a local CellData object
    CellData2P2C()
    {
        for (int i = 0; i < numPhases;i++)
        {
            mobility_[i] = 0.0;
            dv_[i] = 0.0;
        }
        volumeError_ = 0.;
        errorCorrection_= 0.;
        dv_dp_ = 0.;
        volumeDerivativesAvailable_ = false;
        wasRefined_ = false;
        globalIdx_ = 0;
        perimeter_ = 0.;
    }
    //! Acess to flux data, representing information living on the intersections
    FluxData& fluxData()
    {
        return fluxData_;
    }
    //! Constant acess to flux data, representing information living on the intersections
    const FluxData& fluxData() const
    {
        return fluxData_;
    }

    /*! \name Acess to primary variables */
    //@{
    /*!
     * \brief Acess to the phase pressure
     * \param phaseIdx index of the Phase
     */
    Scalar pressure(int phaseIdx)
    {
        return fluidState_->pressure(phaseIdx);
    }

    /*!
     * \brief Acess to the phase pressure
     * \param phaseIdx index of the Phase
     */
    const Scalar pressure(int phaseIdx) const
    {
        return fluidState_->pressure(phaseIdx);
    }

    /*!
     * \brief Returns the total mass concentration of a component \f$\mathrm{[kg/m^3]}\f$.
     *
     * This is equivalent to the sum of the component concentrations for all
     * phases multiplied with the phase density.
     *
     * \param compIdx the index of the component
     */
    const Scalar totalConcentration(int compIdx) const
    {
        return massConcentration_[compIdx];
    }

    //! \copydoc CellData2P2C::totalConcentration()
    const Scalar massConcentration(int compIdx) const
    {
        return massConcentration_[compIdx];
    }

    /*!
     * \brief Sets the total mass concentration of a component \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param compIdx index of the Component
     * \param value Value to be stored
     */
    void setTotalConcentration(int compIdx, Scalar value)
    {
        massConcentration_[compIdx] = value;
    }

    //! \copydoc CellData2P2C::setTotalConcentration()
    void setMassConcentration(int compIdx, Scalar value)
    {
        massConcentration_[compIdx] = value;
    }
    /*!
     * \brief Calculate the total mass concentration of a component \f$\mathrm{[kg/m^3]}\f$
     * for a given porosity (within the initialization procedure).
     * \param porosity Porosity
     */
    void calculateMassConcentration(Scalar porosity)
    {
        massConcentration_[wCompIdx] =
                porosity * (massFraction(wPhaseIdx,wCompIdx) * saturation(wPhaseIdx) * density(wPhaseIdx)
              + massFraction(nPhaseIdx,wCompIdx) * saturation(nPhaseIdx) * density(nPhaseIdx));
        massConcentration_[nCompIdx] =
                porosity * (massFraction(wPhaseIdx,nCompIdx) * saturation(wPhaseIdx) * density(wPhaseIdx)
              + massFraction(nPhaseIdx,nCompIdx) * saturation(nPhaseIdx) * density(nPhaseIdx));
    }
    //@}


    /*! \name Acess to secondary variables */
    //@{
    /*!
     * \brief Return phase mobilities
     * @param phaseIdx index of the Phase
     */
    const Scalar& mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Set phase mobilities
     * \param phaseIdx index of the Phase
     * \param value Value to be stored
     */
    void setMobility(int phaseIdx, Scalar value)
    {
        mobility_[phaseIdx]=value;
    }

    /*!
     * \brief Return the volume error [-].
     * This quantity stands for the deviation of real fluid volume
     * to available pore space.
     * \f$ \epsilon = v_{real} - \phi\f$.
     */
    Scalar& volumeError()
    {
        return volumeError_;
    }

    //! Return the volume error [-].
    const Scalar& volumeError() const
    {
        return volumeError_;
    }

    /*!
     * \brief Return the error Correction
     * This quantifies the damped error that actually
     * entered the pressure equation: Damped Error
     * from last time-step times last time step size
     */
    Scalar& errorCorrection()
    {
        return errorCorrection_;
    }

    //! Return the error Correction
    const Scalar& errorCorrection() const
    {
        return errorCorrection_;
    }

    /*!
     * \brief Return the derivative of specific volume w.r.t. pressure
     * For details, see description of FVPressureCompositional<TypeTag>::volumeDerivatives()
     */
    Scalar& dv_dp()
    {
        return dv_dp_;
    }

    //! Return the derivative of specific volume w.r.t. pressure
    const Scalar& dv_dp() const
    {
        return dv_dp_;
    }

    /*!
     * \brief Return the derivative of spec. volume w.r.t. change of mass
     * For details, see description of FVPressureCompositional<TypeTag>::volumeDerivatives()
     * \param compIdx index of the Component
     */
    Scalar& dv(int compIdx)
    {
        return dv_[compIdx];
    }

    //! \copydoc dv()
    const Scalar& dv(int compIdx) const
    {
        return dv_[compIdx];
    }

    /*!
     * \brief Return cell perimeter (as weithing function)
     * The cell perimeter is used in combination with the face Area as a
     * weighting of the volume integral in the pressure equation.
     */
    Scalar& perimeter()
    {
        return perimeter_;
    }
    //! Return cell perimeter (as weithing function)
    const Scalar& perimeter() const
    {
        return perimeter_;
    }

    /*** b) from fluidstate ***/

    //! DOC ME!
    const Scalar saturation(int phaseIdx) const
    {
        return fluidState_->saturation(phaseIdx);
    }

    //! DOC ME!
    const Scalar viscosity(int phaseIdx) const
    {
        return fluidState_->viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    const Scalar capillaryPressure() const
    {
        return fluidState_->pressure(nPhaseIdx) - fluidState_->pressure(wPhaseIdx);
    }

    //! DOC ME!
    const Scalar density(int phaseIdx) const
    {
        return (fluidState_->density(phaseIdx));
    }

    //! DOC ME!
    const Scalar massFraction(int phaseIdx, int compIdx) const
    {
        return fluidState_->massFraction(phaseIdx, compIdx);
    }

    //! DOC ME!
    const Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        return fluidState_->moleFraction(phaseIdx, compIdx);
    }

    //! DOC ME!
    const Scalar temperature(int phaseIdx) const
    {
        return fluidState_->temperature(phaseIdx);
    }

    //! DOC ME!
    const Scalar phaseMassFraction(int phaseIdx) const
    {
        return fluidState_->phaseMassFraction(phaseIdx);
    }
    //@}

    //! Returns a reference to the cells fluid state
    const FluidState& fluidState() const
    {
        return *fluidState_;
    }

    /*!
     * \brief Allows manipulation of the cells fluid state
     * Fluidstate is stored as a pointer, initialized as a null-pointer.
     * Ensure that if no FluidState is present, a new one is created.
     */
    FluidState& manipulateFluidState()
    {
        if(!fluidState_)
            fluidState_ = std::make_shared<FluidState>();
        return *fluidState_;
    }

    //! stores this cell datas index, only for debugging purposes!!
    int& globalIdx()
    { return globalIdx_;}
    //! Indicates if volume derivatives are computed and available
    bool hasVolumeDerivatives() const
    { return volumeDerivativesAvailable_;}
    //! Specifies that volume derivatives are computed and available
    void confirmVolumeDerivatives()
    { volumeDerivativesAvailable_ = true;}
    //! Specifies if volume derivatives are computed and available
    void volumeDerivativesAvailable(bool value)
    { volumeDerivativesAvailable_ = value;}
    //! Resets the cell data after a timestep was completed: No volume derivatives yet available
    void reset()
    {
        wasRefined_ = false;
        volumeDerivativesAvailable_ = false;
    }
    //! Indicates if current cell was refined at this time step
    bool& wasRefined()
    {   return wasRefined_;}
    //!\copydoc wasRefined()
    const bool& wasRefined() const
    {   return wasRefined_;}

    /*!
     * \brief Indicates if current cell is the upwind cell for a given interface
     * \param indexInInside Local face index seen from current cell
     * \param phaseIdx The index of the phase
     */
    const bool& isUpwindCell(int indexInInside, int phaseIdx) const
    {
        return fluxData_.isUpwindCell(indexInInside, phaseIdx);
    }
    /*!
     * \brief Specifies if current cell is the upwind cell for a given interface
     * \param indexInInside Local face index seen from current cell
     * \param phaseIdx The index of the phase
     * \param value Value: true (->outflow) or false (-> inflow)
     */
    void setUpwindCell(int indexInInside, int phaseIdx, bool value)
    {
        fluxData_.setUpwindCell(indexInInside, phaseIdx, value);
    }

};
}
#endif
