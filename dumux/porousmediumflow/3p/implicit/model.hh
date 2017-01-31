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
 * \brief Adaption of the fully implicit scheme to the three-phase flow model.
 *
 * The model is designed for simulating three fluid phases with water, gas, and
 * a liquid contaminant (NAPL - non-aqueous phase liquid)
 */
#ifndef DUMUX_3P_MODEL_HH
#define DUMUX_3P_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePModel
 * \brief Adaption of the fully implicit scheme to the three-phase flow model.
 *
 * This model implements three-phase flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$
 * The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum.
 *
 * By inserting this into the equations for the conservation of the
 * components, the well-known multiphase flow equation is obtained.
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are gas phase pressure \f$p_g\f$,
 * water saturation \f$S_w\f$ and NAPL saturation \f$S_n\f$.
 */
template<class TypeTag>
class ThreePModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);;

    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
    };

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
        vtkOutputModule.addSecondaryVariable("sw", [](const VolumeVariables& v){ return v.saturation(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("sn", [](const VolumeVariables& v){ return v.saturation(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("sg", [](const VolumeVariables& v){ return v.saturation(gPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pn", [](const VolumeVariables& v){ return v.pressure(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pg", [](const VolumeVariables& v){ return v.pressure(gPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("rhow", [](const VolumeVariables& v){ return v.density(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("rhon", [](const VolumeVariables& v){ return v.density(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("rhog", [](const VolumeVariables& v){ return v.density(gPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
        vtkOutputModule.addSecondaryVariable("permeability", [](const VolumeVariables& v){ return v.permeability(); });
        vtkOutputModule.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });
//         vtkOutputModule.addSecondaryVariable("mobW", [](const VolumeVariables& v){ return v.mobility(wPhaseIdx); });
//         vtkOutputModule.addSecondaryVariable("mobN", [](const VolumeVariables& v){ return v.mobility(nPhaseIdx); });
//         vtkOutputModule.addSecondaryVariable("mobG", [](const VolumeVariables& v){ return v.mobility(gPhaseIdx); });
        // TODO: these lines are commented-out in order to comply with the "old" reference solution;
        // can be changed some time, as well as the parameter names
    }
};

}

#include "propertydefaults.hh"

#endif
