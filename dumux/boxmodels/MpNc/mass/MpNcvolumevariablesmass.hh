/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
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
 * \brief Contains the mass conservation part of the volume variables
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_MASS_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_MASS_HH

#include <dumux/material/fluidstates/genericfluidstate.hh>
#include <dumux/material/fluidstates/equilibriumfluidstate.hh>

namespace Dumux
{
/*!
 * \brief The compositional part of the volume variables if chemical
 *        equilibrium _is_ assumed
 */
template <class TypeTag, bool enableKinetic /* = false */>
class MPNCVolumeVariablesMass
{
    static_assert(!enableKinetic,
                  "No kinetic mass transfer module included, "
                  "but kinetic mass transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CompositionFromFugacitiesSolver)) CompositionFromFugacitiesSolver;
    typedef typename GridView::template Codim<0>::Entity Element;


    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { fug0Idx = Indices::fug0Idx };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:
    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state.
     *
     * If chemical equilibrium is assumed, we use the fluid state
     * which saves some memory.
     */
    typedef EquilibriumFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update composition of all phases in the mutable
     *        parameters from the primary variables.
     */
    template <class MutableParams>
    void update(MutableParams &mutParams,
                const PrimaryVariables &priVars,
                const VolumeVariables *hint,

                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx)
    {
        ComponentVector fug;
        // retrieve component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            fug[compIdx] = priVars[fug0Idx + compIdx];

        // calculate phase compositions
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // initial guess
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar x_ij = 1.0/numComponents;
                if (hint)
                    // use the hint for the initial mole fraction!
                    x_ij = hint->fluidState().moleFrac(phaseIdx, compIdx);

                // set initial guess of the component's mole fraction
                mutParams.setMoleFrac(phaseIdx,
                                      compIdx,
                                      x_ij);

            }

            // calculate the phase composition from the component
            // fugacities
            if (!hint)
                // if no hint was given, we ask the constraint solver
                // to make an initial guess
                CompositionFromFugacitiesSolver::guessInitial(mutParams, phaseIdx, fug);
            CompositionFromFugacitiesSolver::solve(mutParams, phaseIdx, fug);

            /*
              std::cout << "final composition: " << FluidSystem::phaseName(phaseIdx) << "="
              << mutParams.moleFrac(phaseIdx, 0) << " "
              << mutParams.moleFrac(phaseIdx, 1) << " "
              << mutParams.moleFrac(phaseIdx, 2) << " "
              << mutParams.moleFrac(phaseIdx, 3) << " "
              << mutParams.moleFrac(phaseIdx, 4) << " "
              << mutParams.moleFrac(phaseIdx, 5) << " "
              << mutParams.moleFrac(phaseIdx, 6) << "\n";
            */

        }
    }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
    }
};

} // end namepace

#endif
