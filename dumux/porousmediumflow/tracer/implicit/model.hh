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
 * \brief Base class for all models which use the fully implicit tracer model.
 *        Adaption of the fully implicit scheme to the tracer transport model.
 */

#ifndef DUMUX_TRACER_MODEL_HH
#define DUMUX_TRACER_MODEL_HH

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TracerModel
 * \brief Adaption of the fully implicit scheme to the tracer transport model.
 *
 * This model implements a transport of a tracer, where density and viscosity of the
 * fluid phase in which the tracer gets transported are not affected by the tracer.
 * The velocity field is a given spatial parameter.
 * The model is mainly useful for fast computations on given or precomputed
 * velocity fields and thus generally makes sense only in combination with an incompressible
 * one-phase flow velocity field or analytically given / artificial fields. However, reactions
 * between multiple tracer components can be implemented.
 *
 * The transport of the components \f$\kappa \in \{ a, b, c, ... \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa {\textbf v_f}
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme (TPFA or MPFA) as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables the mole or mass fraction of dissolved components \f$x\f$.
 * Note that the tracer model is always considered non-isothermal.
 * The velocity output is fully compatible with the tracer model if you want to write the velocity field to vtk.
 */

template<class TypeTag >
class TracerModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

public:

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem& problem)
    {
        ParentType::init(problem);

        // register standardized vtk output fields
        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addSecondaryVariable("x_" + std::string(FluidSystem::componentName(0)),
                                             [](const VolumeVariables& v){ return v.moleFraction(0, 0); });
        vtkOutputModule.addSecondaryVariable("X_" + std::string(FluidSystem::componentName(0)),
                                             [](const VolumeVariables& v){ return v.massFraction(0, 0); });
        vtkOutputModule.addSecondaryVariable("rho", [](const VolumeVariables& v){ return v.density(); });
    }
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
