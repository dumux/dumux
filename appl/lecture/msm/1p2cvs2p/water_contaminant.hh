/*****************************************************************************
 *   Copyright (C) 2010 by Klaus Mosthaf                                    *
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
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
/*!
 * \file
 *
 * \brief A fluid system with one phase and an arbitrary number of components.
 */
#ifndef WATERCONTAMINANT_HH
#define WATERCONTAMINANT_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/boxmodels/1p2c/1p2cproperties.hh>

namespace Dumux
{

namespace Properties
{
    NEW_PROP_TAG(Scalar);
    NEW_PROP_TAG(OnePTwoCIndices);
};

/*!
 * \brief A fluid system with one phase and an arbitrary number of components.
 */
template <class TypeTag, bool verbose=true>
class WaterContaminant
{
    typedef WaterContaminant<TypeTag, verbose> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

public:
    enum {
        // component indices
        waterIdx = 0,
        contaminantIdx = 1
    };

    static void init()
    {}

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        switch(compIdx)
        {
        case waterIdx:
            return "Water";
        case contaminantIdx:
            return "Contaminant";
        default:
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    /*!
     * \brief Return the molar mass of a component [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    { return 18e-3; }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 1000;//constDensity_; // in [kg /m^3]

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the dynamic viscosity of a phase.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        if (phaseIdx == 0)
            return 1e-03;

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficent of a component in a phase.
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
    {
    //TODO: return diffCoefficient_;
        return  1.0e-9; // in [m^2/s]
    //TODO: The example is very bad! The numerical diffusion is very high, so that the diffusion/dispersion coefficients nearly do not have any influence!
    }


/*    WaterContaminant( )
    {
        //load interface-file
        Dumux::InterfaceFluidProperties interfaceFluidProps("interface1p2c.xml");

        diffCoefficient_ = interfaceFluidProps.IFP_MolecularDiffusionCoefficient;
    }*/


//private:
    static Scalar diffCoefficient_;
};

} // end namepace

#endif
