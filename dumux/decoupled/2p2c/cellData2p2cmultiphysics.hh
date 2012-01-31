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
#ifndef DUMUX_ELEMENTDATA2P2C_MULTYPHYSICS_HH
#define DUMUX_ELEMENTDATA2P2C_MULTYPHYSICS_HH

#include "2p2cproperties.hh"
#include "cellData2p2c.hh"
#include <dumux/decoupled/2p2c/pseudo1p2cfluidstate.hh>

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
class CellData2P2Cmultiphysics : public CellData2P2C<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
//    typedef FluxData2P2C<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;
    typedef PseudoOnePTwoCFluidState<TypeTag> SimpleFluidState;

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
    enum
    {
        complex = 0, simple = 1
    };
private:
    int subdomain_;
    int fluidStateType_;
    SimpleFluidState* simpleFluidState_;

//    FluxData fluxData_;
public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    CellData2P2Cmultiphysics() : CellData2P2C<TypeTag>(),
        subdomain_(2), fluidStateType_(complex), simpleFluidState_(0)
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

    void setSimpleFluidState(SimpleFluidState& simpleFluidState)
    {
        assert (this->subdomain() != 2);
        fluidStateType_ = simple;
        simpleFluidState_ = &simpleFluidState;
    }

    SimpleFluidState& manipulateSimpleFluidState()
    {
        assert (this->subdomain() != 2);
        fluidStateType_ = simple;
        if(!simpleFluidState_)
            simpleFluidState_ = new SimpleFluidState;
        return *simpleFluidState_;
    }

//    const FluidState& fluidState() const
//    {
//        return *fluidState_;
//    }

    FluidState& manipulateFluidState()
    {
        assert(this->subdomain() == 2);
        fluidStateType_ = complex;
        if(!this->fluidState_)
            this->fluidState_ = new FluidState;
        return *this->fluidState_;
    }

    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////
    Scalar pressure(int phaseIdx)
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->pressure(phaseIdx);
        }
        else
            return this->fluidState_->pressure(phaseIdx);
    }

    const Scalar pressure(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->pressure(phaseIdx);
        }
        else
            return this->fluidState_->pressure(phaseIdx);
    }

    //! Return saturation vector
    void setTotalConcentration(int compIdx, Scalar value)
    {
        if(fluidStateType_ == simple)
            return simpleFluidState_->setMassConcentration(compIdx, value);
        else
            return this->fluidState_->setMassConcentration(compIdx, value);
    }

    const Scalar totalConcentration(int compIdx) const
    {
        if(fluidStateType_ == simple)
            return simpleFluidState_->massConcentration(compIdx);
        else
            return this->fluidState_->massConcentration(compIdx);
    }
    const Scalar massConcentration(int compIdx) const
    {
        if(fluidStateType_ == simple)
            return simpleFluidState_->massConcentration(compIdx);
        else
            return this->fluidState_->massConcentration(compIdx);
    }

    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////

    //! Return subdomain information
    int& subdomain()
    {
        return subdomain_;
    }
    const int& subdomain() const
    {
        return subdomain_;
    }

    /*** b) from fluidstate ***/

    //! Return saturation vector
    void setSaturation(int phaseIdx, Scalar value)
    {
        assert(this->subdomain() == 2);
        this->fluidState_->setSaturation(phaseIdx, value);
    }

    const Scalar saturation(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->saturation(phaseIdx);
        }
        else
            return this->fluidState_->saturation(phaseIdx);
    }

    //! Return density vector
    void setViscosity(int phaseIdx, Scalar value)
    {
        if(fluidStateType_ == simple)
        {
            assert(phaseIdx == simpleFluidState_->presentPhaseIdx());
            simpleFluidState_->setViscosity(phaseIdx, value);
        }
        else
            this->fluidState_->setViscosity(phaseIdx, value);
    }

    const Scalar viscosity(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            if(phaseIdx != simpleFluidState_->presentPhaseIdx())
                return 0.; // This should only happend for output
            return simpleFluidState_->viscosity(phaseIdx);
        }
        else
            return this->fluidState_->viscosity(phaseIdx);
    }


    //! Return capillary pressure vector
    const Scalar capillaryPressure() const
    {
        if(fluidStateType_ == simple)
            return 0.;
        else
            return this->fluidState_->pressure(nPhaseIdx) - this->fluidState_->pressure(wPhaseIdx);
    }

    //! Return density vector
    const Scalar density(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->density(phaseIdx);
        }
        else
            return this->fluidState_->density(phaseIdx);
    }

    //! Return the mass fraction
    const Scalar massFraction(int phaseIdx, int compIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->massFraction(phaseIdx, compIdx);
        }
        else
            return this->fluidState_->massFraction(phaseIdx, compIdx);
    }

    //! Return the mole fraction
    const Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->moleFraction(phaseIdx, compIdx);
        }
        else
            return this->fluidState_->moleFraction(phaseIdx, compIdx);
    }
    //! Return temperature
    const Scalar temperature(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->temperature(phaseIdx);
        }
        else
            return this->fluidState_->temperature(phaseIdx);
    }

    //! Return phase phase mass fraction
    const Scalar phaseMassFraction(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            if(phaseIdx == simpleFluidState_->presentPhaseIdx())
                return 1.; // phase is present => nu = 1
            else
                return 0.;
        }
        else
            return this->fluidState_->phaseMassFraction(phaseIdx);
    }
};
}
#endif
