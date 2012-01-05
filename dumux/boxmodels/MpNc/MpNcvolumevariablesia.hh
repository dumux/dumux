// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Philipp Nuske                                *
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_IA_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_IA_HH

#include <dumux/boxmodels/MpNc/MpNcproperties.hh>

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
template <class TypeTag, bool enableKinetic /* = false */, bool enableKineticEnergy /* = false */>
class MPNCVolumeVariablesIA
{
    static_assert(not enableKinetic and not enableKineticEnergy,
                  "The kinetic energy modules need specific interfacial area "
                  "but no suitable specialization of the IA volume variables module "
                  "has been included.");

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    template <class FluidState, class ParameterCache>
    void update(const VolumeVariables & volVars,
                const FluidState &fluidState,
                const ParameterCache &paramCache,
                const PrimaryVariables &priVars,
                const Problem &problem,
                const Element & element,
                const FVElementGeometry & elemGeom,
                const int scvIdx)
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
