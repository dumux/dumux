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
 * \brief Storage container for discretized data for multi-physics models.
 */
#ifndef DUMUX_ELEMENTDATA2P2C_MULTYPHYSICS_HH
#define DUMUX_ELEMENTDATA2P2C_MULTYPHYSICS_HH

#include <dumux/material/fluidstates/pseudo1p2c.hh>
#include "properties.hh"
#include "celldata.hh"

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Storage container for discretized data for multiphysics models.
 *
 * For multi-physics models, we divide the model in separate sub-domains. Being a cell-based
 * information, this is also stored in the cellData. In addition, a simpler version of a
 * fluidState can be stored in cells being in the simpler subdomain.
 * Hence, acess functions either direct to the full fluidstate, or to the simple fluidstate.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CellData2P2CMultiPhysics : public CellData2P2C<TypeTag>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SimpleFluidState = PseudoOnePTwoCFluidState<Scalar, FluidSystem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        complex = 0, simple = 1
    };
private:
    int subdomain_;
    int fluidStateType_;
    std::shared_ptr<SimpleFluidState> simpleFluidState_;

//    FluxData fluxData_;
public:

    //! Constructor for a local storage object
    CellData2P2CMultiPhysics() : CellData2P2C<TypeTag>(),
        subdomain_(2), fluidStateType_(complex)
    {
    }

    /*! \name Acess to primary variables */
    //@{
    //! DOC ME!
    Scalar pressure(int phaseIdx)
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->pressure(phaseIdx);
        }
        else
            return this->fluidState_->pressure(phaseIdx);
    }
    //! DOC ME!
    const Scalar pressure(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->pressure(phaseIdx);
        }
        else
            return this->fluidState_->pressure(phaseIdx);
    }
    //! DOC ME!
    void setPressure(int phaseIdx, Scalar value)
    {
        if(fluidStateType_ == simple)
            manipulateSimpleFluidState().setPressure(phaseIdx, value);
        else
            manipulateFluidState().setPressure(phaseIdx, value);
    }
    //@}


    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////

    /*!
     * \brief Return subdomain information
     *
     * Acess function to store subdomain information
     */
    int& subdomain()
    {
        return subdomain_;
    }

    /*!
     * \brief Return subdomain information
     *
     * Acess function to get subdomain information
     */
    const int& subdomain() const
    {
        return subdomain_;
    }

    /*!
     * \brief Specify subdomain information and fluidStateType
     */
    void setSubdomainAndFluidStateType(int index)
    {
        subdomain_ = index;
        if(index == 2)
            fluidStateType_ = complex;
        else
            fluidStateType_ = simple;
    }

    /*! \name Acess to secondary variables */
    //@{
    //! DOC ME!
    void setSaturation(int phaseIdx, Scalar value)
    {
        if(fluidStateType_ == simple)
        {
            // saturation is triggered by presentPhaseIdx
            manipulateSimpleFluidState().setPresentPhaseIdx((value==0.) ? nPhaseIdx : wPhaseIdx);
        }
        else
        {
            manipulateFluidState().setSaturation(phaseIdx, value);
            manipulateFluidState().setSaturation(1-phaseIdx, 1.0-value);
        }
    }
    //! DOC ME!
    const Scalar saturation(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->saturation(phaseIdx);
        }
        else
            return this->fluidState_->saturation(phaseIdx);
    }

    //! DOC ME!
    void setViscosity(int phaseIdx, Scalar value)
    {
        if(fluidStateType_ == simple)
        {
            assert(phaseIdx == simpleFluidState_->presentPhaseIdx());
            manipulateSimpleFluidState().setViscosity(phaseIdx, value);
        }
        else
            manipulateFluidState().setViscosity(phaseIdx, value);
    }
    //! DOC ME!
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


    //! DOC ME!
    const Scalar capillaryPressure() const
    {
        if(fluidStateType_ == simple)
            return simpleFluidState_->pressure(nPhaseIdx) - simpleFluidState_->pressure(wPhaseIdx);
        else
            return this->fluidState_->pressure(nPhaseIdx) - this->fluidState_->pressure(wPhaseIdx);
    }

    //! DOC ME!
    const Scalar density(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->density(phaseIdx);
        }
        else
            return this->fluidState_->density(phaseIdx);
    }

    //! DOC ME!
    const Scalar massFraction(int phaseIdx, int compIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->massFraction(phaseIdx, compIdx);
        }
        else
            return this->fluidState_->massFraction(phaseIdx, compIdx);
    }

    //! DOC ME!
    const Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->moleFraction(phaseIdx, compIdx);
        }
        else
            return this->fluidState_->moleFraction(phaseIdx, compIdx);
    }
    //! DOC ME!
    const Scalar temperature(int phaseIdx) const
    {
        if(fluidStateType_ == simple)
        {
            return simpleFluidState_->temperature(phaseIdx);
        }
        else
            return this->fluidState_->temperature(phaseIdx);
    }

    //! DOC ME!
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

    //! Returns a reference to the cells simple fluid state
    const SimpleFluidState& simpleFluidState() const
    {
        assert(simpleFluidState_);
        return *simpleFluidState_;
    }

    //@}

    /*!
     * \brief Set a simple fluidstate for a cell in the simple domain
     * Uses a simplified fluidstate with less storage capacity
     * and functionality.
     * Makes shure the fluidStateType_ flag is set appropriately in this
     * cell.
     * \param simpleFluidState A fluidstate storing a 1p2c mixture
     */
    void setSimpleFluidState(SimpleFluidState& simpleFluidState)
    {
        assert (this->subdomain() != 2);
        fluidStateType_ = simple;
        simpleFluidState_ = &simpleFluidState;
    }
    //! Manipulates a simple fluidstate for a cell in the simple domain
    SimpleFluidState& manipulateSimpleFluidState()
    {
        fluidStateType_ = simple;
        if(this->fluidState_)
        {
            this->fluidState_.template reset<FluidState>(0);
        }

        if(!simpleFluidState_)
            simpleFluidState_ = std::make_shared<SimpleFluidState>();
        return *simpleFluidState_;
    }

    /*!
     * \brief Allows manipulation of the complex fluid state
     *
     * Fluidstate is stored as a pointer, initialized as a null-pointer.
     * Enshure that if no FluidState is present, a new one is created.
     * Also enshure that we are in the complex subdomain.
     */
    FluidState& manipulateFluidState()
    {
        fluidStateType_ = complex;
        if(simpleFluidState_)
        {
            simpleFluidState_.template reset<SimpleFluidState>(0);
        }

        if(!this->fluidState_)
        {
            this->fluidState_ = std::make_shared<FluidState>();
            // properly initialize pressure, since it is evaluated later:
            this->fluidState_->setPressure(wPhaseIdx, 1e5);
            this->fluidState_->setPressure(nPhaseIdx, 1e5);
        }

        return *this->fluidState_;
    }

    //! Returns the type of the fluidState
    bool fluidStateType() const
    { return fluidStateType_;}

};

} // end namespace Dumux
#endif
