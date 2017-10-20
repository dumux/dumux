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
* \brief  Base class for all models which use the single-phase, n-component fully implicit model.
*         Adaption of the fully implicit model to the one-phase n-component flow model.
*/

#ifndef DUMUX_1PNC_MODEL_HH
#define DUMUX_1PNC_MODEL_HH

#include <dumux/porousmediumflow/nonisothermal/implicit/model.hh>

#include "properties.hh"


namespace Dumux
{
/*!
 * \ingroup OnePNCModel
 * \brief Adaption of the fully implicit scheme to the
 *        one-phase n-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists
 * of n components, using a standard Darcy approach as the equation for the
 * conservation of momentum:
 \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \phi\frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a, ... \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p -
 \varrho {\textbf g} \right)
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p\f$ and the mole fraction of dissolved components \f$x^\kappa\f$.
 */

template<class TypeTag>
class OnePNCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using NonIsothermalModel = Dumux::NonIsothermalModel<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const int phaseIdx = Indices::phaseIdx;

    enum { dim = GridView::dimension };
    enum {  numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

public:
    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        // register standardized vtk output fields
        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addSecondaryVariable("pressure", [](const VolumeVariables& v){ return v.pressure(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("rho", [](const VolumeVariables& v){ return v.density(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });

        for (int i = 0; i < numComponents; ++i)
           vtkOutputModule.addSecondaryVariable("x_" + std::string(FluidSystem::componentName(i)),
                                                [i](const VolumeVariables& v){ return v.moleFraction(phaseIdx, i); });

        for (int i = 0; i < numComponents; ++i)
           vtkOutputModule.addSecondaryVariable("X_" + std::string(FluidSystem::componentName(i)),
                                                 [i](const VolumeVariables& v){ return v.massFraction(phaseIdx,i); });

        NonIsothermalModel::maybeAddTemperature(vtkOutputModule);
    }
};

}

#include "propertydefaults.hh"

#endif
