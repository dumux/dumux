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
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_IA_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_IA_HH

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>

namespace Dumux
{

/*!
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 *
 * This is the specialization for the cases which do _not_ require
 * specific interfacial area.
 */
template <class TypeTag, bool enableKinetic, int numEnergyEquations>
class MPNCVolumeVariablesIA
{
    static_assert(((numEnergyEquations < 1) && !enableKinetic),
                  "The kinetic energy modules need specific interfacial area "
                  "but no suitable specialization of the IA volume variables module "
                  "has been included.");

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum {dimWorld=GridView::dimensionworld};
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     *
     *      \param volVars The volume variables
     *      \param fluidState Container for all the secondary variables concerning the fluids
     *      \param paramCache Container for cache parameters
     *      \param priVars The primary Variables
     *      \param problem The problem
     *      \param element The finite element
     *      \param fvGeometry The finite-volume geometry in the fully implicit scheme
     *      \param scvIdx The index of the sub-control volumete element
     *
     */
    void update(const VolumeVariables & volVars,
                const FluidState &fluidState,
                const ParameterCache &paramCache,
                const PrimaryVariables &priVars,
                const Problem &problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx)
    {
    }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    { }
};

} // namespace Dumux

#endif
