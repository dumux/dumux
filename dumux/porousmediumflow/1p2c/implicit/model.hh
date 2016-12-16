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
 * \brief Base class for all models which use the single-phase,
 *        two-component fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase two-component flow model.
 */

#ifndef DUMUX_ONEP_TWOC_MODEL_HH
#define DUMUX_ONEP_TWOC_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCModel
 * \brief Adaption of the fully implicit scheme to the one-phase two-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists of two components,
 * using a standard Darcy
 * approach as the equation for the conservation of momentum:
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
 * The transport of the components \f$\kappa \in \{ w, a \}\f$ is described by the following equation:
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
 * The primary variables are the pressure \f$p\f$ and the mole or mass fraction of dissolved component \f$x\f$.
 */

template<class TypeTag >
class OnePTwoCModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int phaseIdx = Indices::phaseIdx;
    static const bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);

public:
    /*!
     * \brief \copybrief ImplicitModel::addOutputVtkFields
     *
     * Specialization for the OnePTwoCModel, adding pressure,
     * mass and mole fractions, and the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        // create the required scalar fields
        unsigned numDofs = this->numDofs();
        auto& pressure = *writer.allocateManagedBuffer(numDofs);
        auto& delp = *writer.allocateManagedBuffer(numDofs);
        auto& moleFraction0 = *writer.allocateManagedBuffer(numDofs);
        auto& moleFraction1 = *writer.allocateManagedBuffer(numDofs);
        auto& massFraction0 = *writer.allocateManagedBuffer(numDofs);
        auto& massFraction1 = *writer.allocateManagedBuffer(numDofs);
        auto& rho = *writer.allocateManagedBuffer(numDofs);
        auto& mu = *writer.allocateManagedBuffer(numDofs);

        auto& velocity = *(writer.template allocateManagedBuffer<double, dimWorld>(numDofs));
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput())
            velocity = 0.0;

        unsigned numElements = this->gridView_().size(0);
        auto& rank = *writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
           int eIdx = this->problem_().model().elementMapper().index(element);

            rank[eIdx] = this->gridView_().comm().rank();

            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bind(element);

            auto elemVolVars = localView(this->curGlobalVolVars());
            elemVolVars.bind(element, fvGeometry, this->curSol());

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                const auto dofIdxGlobal = scv.dofIndex();

                pressure[dofIdxGlobal] = volVars.pressure(phaseIdx);
                delp[dofIdxGlobal] = volVars.pressure(phaseIdx) - 1e5;
                moleFraction0[dofIdxGlobal] = volVars.moleFraction(phaseIdx, 0);
                moleFraction1[dofIdxGlobal] = volVars.moleFraction(phaseIdx, 1);
                massFraction0[dofIdxGlobal] = volVars.massFraction(phaseIdx, 0);
                massFraction1[dofIdxGlobal] = volVars.massFraction(phaseIdx, 1);
                rho[dofIdxGlobal] = volVars.density(phaseIdx);
                mu[dofIdxGlobal] = volVars.viscosity(phaseIdx);
            }

            velocityOutput.calculateVelocity(velocity, elemVolVars, fvGeometry, element, phaseIdx);
        }

        writer.attachDofData(pressure, "P", isBox);
        writer.attachDofData(delp, "delp", isBox);
        if (velocityOutput.enableOutput())
            writer.attachDofData(velocity,  "velocity", isBox, dim);

        writer.attachDofData(moleFraction0, "x_" + std::string(FluidSystem::componentName(0)), isBox);
        writer.attachDofData(moleFraction1, "x_" + std::string(FluidSystem::componentName(1)), isBox);
        writer.attachDofData(massFraction0, "X_" + std::string(FluidSystem::componentName(0)), isBox);
        writer.attachDofData(massFraction1, "X_" + std::string(FluidSystem::componentName(1)), isBox);

        writer.attachDofData(rho, "rho", isBox);
        writer.attachDofData(mu, "mu", isBox);
        writer.attachCellData(rank, "process rank");
    }
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
