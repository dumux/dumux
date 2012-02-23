// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Benjamin Faigle                                   *
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
#ifndef DUMUX_ELEMENTDATA2P2C_HH
#define DUMUX_ELEMENTDATA2P2C_HH

#include "2p2cproperties.hh"
#include <dumux/decoupled/2p2c/dec2p2cfluidstate.hh>
//#include "fluxData2p2c.hh"

/**
 * @file
 * @brief  Storage container for discretized data of the constitutive relations for one element
 * @author Benjamin Faigle
 */

namespace Dumux
{
/*!
 * \ingroup multiphysics multiphase
 */
//! Storage container for discretized data of the constitutive relations for one element
/*! This class stores all cell-centered (FV-Scheme) values for decoupled compositional two-phase flow
 * models that are used by both pressure and transport model. All fluid data are already stored in the
 * fluidstate, so the CellData contains the fluidstate object for the current element.
 * At the moment, the compositional model does not use fluxVariables that are stored on the interfaces.
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CellData2P2C
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
//    typedef FluxData2P2C<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };
protected:
    Scalar mobility_[numPhases];

    //numerical quantities
    Scalar numericalDensity_[numPhases];
    Scalar volumeError_;
    Scalar errorCorrection_;
    Scalar dv_dp_;
    Scalar dv_[numComponents];
    bool volumeDerivativesAvailable_;

    int globalIdx_;
    Scalar perimeter_;

    FluidState* fluidState_;
//    FluxData fluxData_;
public:

    //! Constructor for a local CellData object
    CellData2P2C() :
        fluidState_(0)
    {
        for (int i = 0; i < numPhases;i++)
        {
            mobility_[i] = 0.0;
            numericalDensity_[i] = 0.0;
            mobility_[i] = 0.0;
            dv_[i] = 0.0;
        }
        volumeError_ = 0.;
        errorCorrection_= 0.;
        dv_dp_ = 0.;
        volumeDerivativesAvailable_ = false;
        globalIdx_ = 0;
        perimeter_ = 0.;
    }

//    FluxData& setFluxData()
//    {
//        return fluxData_;
//    }
//
//    const FluxData& fluxData() const
//    {
//        return fluxData_;
//    }


    /*! \name Acess to primary variables */
    //@{
    /*! Acess to the phase pressure
     * @param phaseIdx index of the Phase
     */
    Scalar pressure(int phaseIdx)
    {
        return fluidState_->pressure(phaseIdx);
    }
    /*! Acess to the phase pressure
     * @param phaseIdx index of the Phase
     */
    const Scalar pressure(int phaseIdx) const
    {
        return fluidState_->pressure(phaseIdx);
    }
    /*! Modify the phase pressure
     * @param phaseIdx index of the Phase
     * @param value Value to be srored
     */
    void setPressure(int phaseIdx, Scalar value)
    {
        fluidState_->setPressure(phaseIdx, value);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::massConcentration()
    const Scalar totalConcentration(int compIdx) const
    {
        return fluidState_->massConcentration(compIdx);
    }
    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::massConcentration()
    const Scalar massConcentration(int compIdx) const
    {
        return fluidState_->massConcentration(compIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::setMassConcentration()
    void setTotalConcentration(int compIdx, Scalar value)
    {
        fluidState_->setMassConcentration(compIdx, value);
    }
    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::setMassConcentration()
    void setMassConcentration(int compIdx, Scalar value)
    {
        fluidState_->setMassConcentration(compIdx, value);
    }
    //@}


    /*! \name Acess to secondary variables */
    //@{
    //! Return phase mobilities
    /*
     * @param phaseIdx index of the Phase
     */
    const Scalar& mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }
    //! Set phase mobilities
    /*
     * @param phaseIdx index of the Phase
     * @param value Value to be stored
     */
    void setMobility(int phaseIdx, Scalar value)
    {
        mobility_[phaseIdx]=value;
    }

    //! Return numerical density \f$\mathrm{[kg/m^3]}\f$.
    /** Returns the real fluid density if the whole pore volume
     * was filled. This quantity equals the real density scaled with
     * the volume error.
     * @param phaseIdx index of the Phase
     */
    Scalar& numericalDensity(int phaseIdx)
    {
        return numericalDensity_[phaseIdx];
    }
    //! Return numerical density
    const Scalar& numericalDensity(int phaseIdx) const
    {
        return numericalDensity_[phaseIdx];
    }

    //! Return the volume error [-].
    /** This quantity stands for the deviation of real fluid volume
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

    //! Return the error Correction
    /** This quantifies the damped error that actually
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
    //! Return the derivative of specific volume w.r.t. pressure
    /**
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

    //! Return the derivative of spec. volume w.r.t. change of mass
    /**
     * For details, see description of FVPressureCompositional<TypeTag>::volumeDerivatives()
     * @param compIdx index of the Component
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

    //! Return cell perimeter (as weithing function)
    /*
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

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::setSaturation()
    void setSaturation(int phaseIdx, Scalar value)
    {
        fluidState_->setSaturation(phaseIdx, value);
    }
    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::saturation()
    const Scalar saturation(int phaseIdx) const
    {
        return fluidState_->saturation(phaseIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::setViscosity()
    void setViscosity(int phaseIdx, Scalar value)
    {
        fluidState_->setViscosity(phaseIdx, value);
    }
    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::viscosity()
    const Scalar viscosity(int phaseIdx) const
    {
        return fluidState_->viscosity(phaseIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::capillaryPressure()
    const Scalar capillaryPressure() const
    {
        return fluidState_->pressure(nPhaseIdx) - fluidState_->pressure(wPhaseIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::density()
    const Scalar density(int phaseIdx) const
    {
        return (fluidState_->density(phaseIdx));
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::massFraction()
    const Scalar massFraction(int phaseIdx, int compIdx) const
    {
        return fluidState_->massFraction(phaseIdx, compIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::moleFraction()
    const Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        return fluidState_->moleFraction(phaseIdx, compIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::temperature()
    const Scalar temperature(int phaseIdx) const
    {
        return fluidState_->temperature(phaseIdx);
    }

    //! \copydoc Dumux::DecoupledTwoPTwoCFluidState::phaseMassFraction()
    const Scalar phaseMassFraction(int phaseIdx) const
    {
        return fluidState_->phaseMassFraction(phaseIdx);
    }
    //@}


//    void setFluidState(FluidState& fluidState)
//    DUNE_DEPRECATED
//    {
//        fluidState_ = &fluidState;
//    }

    //! Returns a reference to the cells fluid state
    const FluidState& fluidState() const
    {
        return *fluidState_;
    }

    //! Allows manipulation of the cells fluid state
    /** Fluidstate is stored as a pointer, initialized as a null-pointer.
     * Enshure that if no FluidState is present, a new one is created.
     */
    FluidState& manipulateFluidState()
    {
        if(!fluidState_)
            fluidState_ = new FluidState;
        return *fluidState_;
    }
    //! stores this cell datas index, only for debugging purposes!!
    int& globalIdx()
    { return globalIdx_;}
    //! Indicates if volume derivatives are computed and available
    const bool hasVolumeDerivatives() const
    { return volumeDerivativesAvailable_;}
    //! Specifies that volume derivatives are computed and available
    void confirmVolumeDerivatives()
    { volumeDerivativesAvailable_ = true;}
    //! Resets the cell data after a timestep was completed: No volume derivatives yet available
    void reset()
    {
        volumeDerivativesAvailable_ = false;
    }


};
}
#endif
