// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Contains the mass conservation part of the volume variables.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_MASS_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_MASS_HH

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>

#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux
{
/*!
 * \brief The compositional part of the volume variables if chemical
 *        equilibrium _is_ assumed.
 */
template <class TypeTag, bool enableKinetic /* = false */>
class MPNCVolumeVariablesMass
{
    static_assert(!enableKinetic,
                  "No kinetic mass transfer module included, "
                  "but kinetic mass transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintSolver) ConstraintSolver;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { fug0Idx = Indices::fug0Idx };
    enum { dimWorld = GridView::dimensionworld};
    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:
    /*!
     * \brief Update composition of all phases in the mutable
     *        parameters from the primary variables.
     *
     *        \param fs Container for all the secondary variables concerning the fluids
     *        \param paramCache Container for cache parameters
     *        \param priVars The primary Variables
     *        \param *hint the volume variables, usable for initial guess of composition
     *        \param problem The problem
     *        \param element The finite element
     *        \param fvGeometry The finite-volume geometry in the fully implicit scheme
     *        \param scvIdx The index of the sub-control volume
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const PrimaryVariables &priVars,
                const VolumeVariables *hint,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const unsigned int scvIdx)
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
                    x_ij = hint->fluidState().moleFraction(phaseIdx, compIdx);

                // set initial guess of the component's mole fraction
                fs.setMoleFraction(phaseIdx,
                                   compIdx,
                                   x_ij);

            }

            // calculate the phase composition from the component
            // fugacities
            if (!hint)
                // if no hint was given, we ask the constraint solver
                // to make an initial guess
                ConstraintSolver::guessInitial(fs, paramCache, phaseIdx, fug);
            ConstraintSolver::solve(fs, paramCache, phaseIdx, fug);
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

} // end namespace

#endif
