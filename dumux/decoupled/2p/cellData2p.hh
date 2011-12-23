/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
#ifndef DUMUX_ELEMENTDATA2P_HH
#define DUMUX_ELEMENTDATA2P_HH

#include "2pproperties.hh"
#include "fluxData2p.hh"

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
template<class TypeTag, bool enableCompressibility>
class CellData2P;



template<class TypeTag>
class CellData2P<TypeTag, false>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef FluxData2P<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };
private:
    Scalar saturation_[numPhases];
    Scalar pressure_[numPhases];
    Scalar capillaryPressure_;
    Scalar mobility_[numPhases];
    Scalar fracFlowFunc_[numPhases];
    Scalar update_;

    FluxData fluxData_;
public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    CellData2P() :
        saturation_({0.0, 0.0}), pressure_({0.0, 0.0}), capillaryPressure_(0), mobility_({0.0, 0.0}), fracFlowFunc_({0.0, 0.0})
    {
    }

    FluxData& fluxData()
    {
        return fluxData_;
    }

    const FluxData& fluxData() const
    {
        return fluxData_;
    }

    FluidState& fluidState()
    {
        DUNE_THROW(Dune::NotImplemented,"fluid states not stored in cell data of incompressible models!");
    }

    const FluidState& fluidState() const
    {
        DUNE_THROW(Dune::NotImplemented,"fluid states not stored in cell data of incompressible models!");
    }

    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////


    Scalar pressure(int phaseIdx)
    {
        return pressure_[phaseIdx];
    }

    void setPressure(int phaseIdx, Scalar press)
    {
        pressure_[phaseIdx] = press;
    }

    Scalar& globalPressure()
    {
        return pressure_[wPhaseIdx];
    }

    const Scalar& globalPressure() const
    {
        return pressure_[wPhaseIdx];
    }

    //! Return saturation vector
    Scalar saturation(int phaseIdx)
    {
        return saturation_[phaseIdx];
    }

    void setSaturation(int phaseIdx, Scalar sat)
    {
        saturation_[phaseIdx] = sat;
    }

    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////

    //! Return phase mobilities
    Scalar& mobility(int phaseIdx)
    {
        return mobility_[phaseIdx];
    }

    const Scalar& mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    //! Return phase fractional flow functions
    Scalar& fracFlowFunc(int phaseIdx)
    {
        return fracFlowFunc_[phaseIdx];
    }

    const Scalar& fracFlowFunc(int phaseIdx) const
    {
        return fracFlowFunc_[phaseIdx];
    }

    //! Return capillary pressure vector
    Scalar capillaryPressure()
    {
        return capillaryPressure_;
    }

    void setCapillaryPressure(Scalar pc)
    {
        capillaryPressure_ = pc;
    }

    void setUpdate(Scalar update)
    {
        update_ = update;
    }

    //! Return density vector
    Scalar& volumeCorrection()
    {
        return update_;
    }

    const Scalar& volumeCorrection() const
    {
        return update_;
    }

    //! Return density vector
    Scalar density(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,"density function not implemented in cell data of incompressible models!");
    }
    Scalar density(int phaseIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,"density function not implemented in cell data of incompressible models!");
    }


    //! Return density vector
    Scalar viscosity(int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented,"viscosity function not implemented in cell data of incompressible models!");
    }

    Scalar viscosity(int phaseIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,"viscosity function not implemented in cell data of incompressible models!");
    }
};

template<class TypeTag>
class CellData2P<TypeTag, true>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef FluxData2P<TypeTag> FluxData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };
private:
    Scalar mobility_[numPhases];
    Scalar fracFlowFunc_[numPhases];
    Scalar update_;

    FluxData fluxData_;
    FluidState fluidState_;
public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    CellData2P() :
        mobility_({0.0, 0.0}), fracFlowFunc_({0.0, 0.0})
    {
    }

    FluxData& fluxData()
    {
        return fluxData_;
    }

    const FluxData& fluxData() const
    {
        return fluxData_;
    }

    FluidState& fluidState()
    {
        return fluidState_;
    }

    const FluidState& fluidState() const
    {
        return fluidState_;
    }

    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////


    Scalar pressure(int phaseIdx)
    {
        return fluidState_.pressure(phaseIdx);
    }

    Scalar pressure(int phaseIdx) const
    {
        return fluidState_.pressure(phaseIdx);
    }

    void setPressure(int phaseIdx, Scalar press)
    {
        DUNE_THROW(Dune::NotImplemented,"setPressure function defined for compressible models!");
    }

    Scalar& globalPressure()
    {
        DUNE_THROW(Dune::NotImplemented,"no global pressure defined for compressible models!");
    }

    const Scalar& globalPressure() const
    {
        DUNE_THROW(Dune::NotImplemented,"no global pressure defined for compressible models!");
    }


    //! Return saturation vector
    Scalar saturation(int phaseIdx)
    {
        return fluidState_.saturation(phaseIdx);
    }

    Scalar saturation(int phaseIdx) const
    {
        return fluidState_.saturation(phaseIdx);
    }

    void setSaturation(int phaseIdx, Scalar sat)
    {
        DUNE_THROW(Dune::NotImplemented,"setSaturation function defined for compressible models!");
    }

    //////////////////////////////////////////////////////////////
    // functions returning the vectors of secondary variables
    //////////////////////////////////////////////////////////////

    //! Return phase mobilities
    Scalar& mobility(int phaseIdx)
    {
        return mobility_[phaseIdx];
    }

    const Scalar& mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    //! Return phase fractional flow functions
    Scalar& fracFlowFunc(int phaseIdx)
    {
        return fracFlowFunc_[phaseIdx];
    }

    const Scalar& fracFlowFunc(int phaseIdx) const
    {
        return fracFlowFunc_[phaseIdx];
    }

    //! Return capillary pressure vector
    Scalar capillaryPressure()
    {
        return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx);
    }

    Scalar capillaryPressure() const
    {
        return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx);
    }

    void setCapillaryPressure(Scalar pc)
    {
        DUNE_THROW(Dune::NotImplemented,"no capillary pressure stored for compressible models!");
    }

    //! Return density vector
    Scalar density(int phaseIdx)
    {
        return fluidState_.density(phaseIdx);
    }
    Scalar density(int phaseIdx) const
    {
        return fluidState_.density(phaseIdx);
    }


    //! Return density vector
    Scalar viscosity(int phaseIdx)
    {
        return fluidState_.viscosity(phaseIdx);
    }

    Scalar viscosity(int phaseIdx) const
    {
        return fluidState_.viscosity(phaseIdx);
    }

    void setUpdate(Scalar update)
    {
        update_ = update;
    }

    //! Return density vector
    Scalar& volumeCorrection()
    {
        return update_;
    }

    const Scalar& volumeCorrection() const
    {
        return update_;
    }
};
}
#endif
