/*****************************************************************************
 *   Copyright (C) 2011 by Benjamin Faigle                                   *
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
#ifndef DUMUX_ELEMENTDATA2P2C_HH
#define DUMUX_ELEMENTDATA2P2C_HH

#include "2p2cproperties.hh"
#include <dumux/decoupled/2p2c/dec2p2cfluidstate.hh>
//#include "fluxData2p2c.hh"

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
/*! TODO: The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * @tparam TypeTag The Type Tag
 1*/
template<class TypeTag>
class CellData2P2C
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
//    typedef FluxData2P2C<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents))
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

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    CellData2P2C() :
        mobility_({0.0, 0.0}),
        numericalDensity_({0.0, 0.0}), volumeError_(0.), errorCorrection_(0.),
        dv_dp_(0.), dv_({0.0, 0.0}), volumeDerivativesAvailable_(false),
        globalIdx_(0), perimeter_(0.),fluidState_(0)
    {
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

    void setFluidState(FluidState& fluidState)
    DUNE_DEPRECATED
    {
        fluidState_ = &fluidState;
    }

    const FluidState& fluidState() const
    {
        return *fluidState_;
    }

    FluidState& manipulateFluidState()
    {
        if(!fluidState_)
            fluidState_ = new FluidState;
        return *fluidState_;
    }

    int& globalIdx()
    { return globalIdx_;}
    const bool hasVolumeDerivatives() const
    { return volumeDerivativesAvailable_;}
    void confirmVolumeDerivatives()
    { volumeDerivativesAvailable_ = true;}
    void reset()
    {
        volumeDerivativesAvailable_ = false;
        //TODO:reset flux stuff!
    }
    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////
    Scalar pressure(int phaseIdx)
    {
        return fluidState_->pressure(phaseIdx);
    }

    const Scalar pressure(int phaseIdx) const
    {
        return fluidState_->pressure(phaseIdx);
    }

    void setPressure(int phaseIdx, Scalar value)
    {
        fluidState_->setPressure(phaseIdx, value);
    }

    //! Return saturation vector
    void setTotalConcentration(int compIdx, Scalar value)
    {
        fluidState_->setMassConcentration(compIdx, value);
    }
    void setMassConcentration(int compIdx, Scalar value)
    {
        fluidState_->setMassConcentration(compIdx, value);
    }

    const Scalar totalConcentration(int compIdx) const
    {
        return fluidState_->massConcentration(compIdx);
    }
    const Scalar massConcentration(int compIdx) const
    {
        return fluidState_->massConcentration(compIdx);
    }

    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////


    //! Return phase mobilities
    void setMobility(int phaseIdx, Scalar value)
    {
        mobility_[phaseIdx]=value;
    }

    const Scalar& mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    //! Return numerical density vector
    Scalar& numericalDensity(int phaseIdx)
    {
        return numericalDensity_[phaseIdx];
    }
    const Scalar& numericalDensity(int phaseIdx) const
    {
        return numericalDensity_[phaseIdx];
    }

    //! Return the volume error
    Scalar& volumeError()
    {
        return volumeError_;
    }
    const Scalar& volumeError() const
    {
        return volumeError_;
    }

    //! Return the error Correction
    Scalar& errorCorrection()
    {
        return errorCorrection_;
    }
    //! Return the error Correction
    const Scalar& errorCorrection() const
    {
        return errorCorrection_;
    }
    //! Return the derivative of spec. volume w.r.t. pressure
    Scalar& dv_dp()
    {
        return dv_dp_;
    }
    const Scalar& dv_dp() const
    {
        return dv_dp_;
    }

    //! Return the derivative of spec. volume w.r.t. mass change
    Scalar& dv(int compIdx)
    {
        return dv_[compIdx];
    }
    const Scalar& dv(int compIdx) const
    {
        return dv_[compIdx];
    }

    //! Return cell perimeter (as weithing function)
    Scalar& perimeter()
    {
        return perimeter_;
    }
    const Scalar& perimeter() const
    {
        return perimeter_;
    }
//
//    //! Return subdomain information
//    int& subdomain()
//    {
//        return subdomain_;
//    }
//    const int& subdomain() const
//    {
//        return subdomain_;
//    }
//
    /*** b) from fluidstate ***/

    //! Return saturation vector
    void setSaturation(int phaseIdx, Scalar value)
    {
        fluidState_->setSaturation(phaseIdx, value);
    }

    const Scalar saturation(int phaseIdx) const
    {
        return fluidState_->saturation(phaseIdx);
    }

    //! Return density vector
    void setViscosity(int phaseIdx, Scalar value)
    {
        fluidState_->setViscosity(phaseIdx, value);
    }

    const Scalar viscosity(int phaseIdx) const
    {
        return fluidState_->viscosity(phaseIdx);
    }


    //! Return capillary pressure vector
    const Scalar capillaryPressure() const
    {
        return fluidState_->pressure(nPhaseIdx) - fluidState_->pressure(wPhaseIdx);
    }

    void setCapillaryPressure(Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented,"no capillary pressure stored for compressible models!");
    }

    //! Return density vector
    const Scalar density(int phaseIdx) const
    {
        return (fluidState_->density(phaseIdx));
    }

    //! Return density vector
    const Scalar massFraction(int phaseIdx, int compIdx) const
    {
        return fluidState_->massFraction(phaseIdx, compIdx);
    }

    //! Return density vector
    const Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        return fluidState_->moleFraction(phaseIdx, compIdx);
    }
    //! Return temperature
    const Scalar temperature(int phaseIdx) const
    {
        return fluidState_->temperature(phaseIdx);
    }

    //! Return phase phase mass fraction
    const Scalar phaseMassFraction(int phaseIdx) const
    {
        return fluidState_->phaseMassFraction(phaseIdx);
    }

};
}
#endif
