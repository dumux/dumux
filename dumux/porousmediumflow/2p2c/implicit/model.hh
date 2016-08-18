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
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase two-component fully implicit model.
 */
#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include "properties.hh"
#include "indices.hh"
#include "primaryvariableswitch.hh"

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase two-component fully implicit model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha \frac{M^\kappa}{M_\alpha} x_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha \frac{M^\kappa}{M_\alpha} x_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 (\textbf{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{ D_{\alpha,\text{pm}}^\kappa \varrho_{\alpha} \frac{M^\kappa}{M_\alpha}
 \textbf{grad} x^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$x^\kappa_w + x^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 * Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^a_w<x^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^w_n<x^w_{n,max})\f$</li>
 * </ul>
 */

template<class TypeTag>
class TwoPTwoCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    // the parent class needs to access the variable switch
    friend typename GET_PROP_TYPE(TypeTag, BaseModel);

    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        pwsn = TwoPTwoCFormulation::pwsn,
        pnsw = TwoPTwoCFormulation::pnsw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    template<typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), void>::type
    newtonEndStep()
    {
        // \todo resize volvars vector if grid was adapted

        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol());

        // update the secondary variables if global caching is enabled
        // \note we only updated if phase presence changed as the volume variables
        //       are already updated once by the switch
        for (const auto& element : elements(this->problem_().gridView()))
        {
            // make sure FVElementGeometry & vol vars are bound to the element
            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            if (switchFlag_)
            {
                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();
                    if (priVarSwitch_().wasSwitched(dofIdxGlobal))
                    {
                        this->nonConstCurGlobalVolVars().volVars(dofIdxGlobal).update(this->curSol()[dofIdxGlobal],
                                                                                      this->problem_(),
                                                                                      element,
                                                                                      scv);
                    }
                }
            }

            // handle the boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip the rest
                if (!scvf.boundary())
                    continue;

                // check if boundary is a pure dirichlet boundary
                const auto bcTypes = this->problem_().boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    const auto insideScvIdx = scvf.insideScvIdx();
                    const auto& insideScv = fvGeometry.scv(insideScvIdx);
                    const auto dirichletPriVars = this->problem_().dirichlet(element, scvf);

                    this->nonConstCurGlobalVolVars().volVars(scvf.outsideScvIdx()).update(dirichletPriVars, this->problem_(), element, insideScv);
                }
            }
        }
    }

    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    template<typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), void>::type
    newtonEndStep()
    {
        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol());
    }

    /*!
     * \brief Called by the update() method if applying the Newton
     *        method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        switchFlag_ = false;
        priVarSwitch_().resetPhasePresence();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();

        // update the phase state
        priVarSwitch_().updateOldPhasePresence();
        switchFlag_ = false;
    }

    /*!
     * \brief Returns true if the primary variables were switched for
     *        at least one dof after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                                MultiWriter &writer)
    {
        using ScalarField = Dune::BlockVector<Dune::FieldVector<double, 1> >;
        // using VectorField = Dune::BlockVector<Dune::FieldVector<double, dim>>;

        // get the number of degrees of freedom
        auto numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *sN    = writer.allocateManagedBuffer(numDofs);
        ScalarField *sW    = writer.allocateManagedBuffer(numDofs);
        ScalarField *pn    = writer.allocateManagedBuffer(numDofs);
        ScalarField *pw    = writer.allocateManagedBuffer(numDofs);
        ScalarField *pc    = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoW  = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoN  = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobW  = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobN = writer.allocateManagedBuffer(numDofs);
        ScalarField *phasePresence = writer.allocateManagedBuffer(numDofs);
        ScalarField *massFrac[numPhases][numComponents];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                massFrac[phaseIdx][compIdx] = writer.allocateManagedBuffer(numDofs);
        ScalarField *moleFrac[numPhases][numComponents];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        moleFrac[phaseIdx][compIdx] = writer.allocateManagedBuffer(numDofs);
        ScalarField *temperature = writer.allocateManagedBuffer(numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        // VectorField *velocityN = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // VectorField *velocityW = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        // if (velocityOutput.enableOutput()) // check if velocity output is demanded
        // {
        //     // initialize velocity fields
        //     for (unsigned int i = 0; i < numDofs; ++i)
        //     {
        //         (*velocityN)[i] = Scalar(0);
        //         (*velocityW)[i] = Scalar(0);
        //     }
        // }

        auto numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            auto eIdxGlobal = this->problem_().elementMapper().index(element);
            (*rank)[eIdxGlobal] = this->gridView_().comm().rank();

            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(this->curGlobalVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->curSol());

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                auto dofIdxGlobal = scv.dofIndex();

                (*sN)[dofIdxGlobal]    = volVars.saturation(nPhaseIdx);
                (*sW)[dofIdxGlobal]    = volVars.saturation(wPhaseIdx);
                (*pn)[dofIdxGlobal]    = volVars.pressure(nPhaseIdx);
                (*pw)[dofIdxGlobal]    = volVars.pressure(wPhaseIdx);
                (*pc)[dofIdxGlobal]    = volVars.capillaryPressure();
                (*rhoW)[dofIdxGlobal]  = volVars.density(wPhaseIdx);
                (*rhoN)[dofIdxGlobal]  = volVars.density(nPhaseIdx);
                (*mobW)[dofIdxGlobal]  = volVars.mobility(wPhaseIdx);
                (*mobN)[dofIdxGlobal]  = volVars.mobility(nPhaseIdx);
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*massFrac[phaseIdx][compIdx])[dofIdxGlobal] = volVars.massFraction(phaseIdx, compIdx);
                        Valgrind::CheckDefined((*massFrac[phaseIdx][compIdx])[dofIdxGlobal][0]);
                    }
                }
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*moleFrac[phaseIdx][compIdx])[dofIdxGlobal] = volVars.moleFraction(phaseIdx, compIdx);
                        Valgrind::CheckDefined((*moleFrac[phaseIdx][compIdx])[dofIdxGlobal][0]);
                    }
                }
                (*poro)[dofIdxGlobal]  = volVars.porosity();
                (*temperature)[dofIdxGlobal] = volVars.temperature();
                (*phasePresence)[dofIdxGlobal]= priVarSwitch().phasePresence(dofIdxGlobal);
            }

            // // velocity output
            // velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
            // velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);

        } // loop over elements

        writer.attachDofData(*sN,     "Sn", isBox);
        writer.attachDofData(*sW,     "Sw", isBox);
        writer.attachDofData(*pn,     "pn", isBox);
        writer.attachDofData(*pw,     "pw", isBox);
        writer.attachDofData(*pc,     "pc", isBox);
        writer.attachDofData(*rhoW,   "rhoW", isBox);
        writer.attachDofData(*rhoN,   "rhoN", isBox);
        writer.attachDofData(*mobW,   "mobW", isBox);
        writer.attachDofData(*mobN,   "mobN", isBox);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "X_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachDofData(*massFrac[phaseIdx][compIdx], oss.str(), isBox);
            }
        }
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "x_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachDofData(*moleFrac[phaseIdx][compIdx], oss.str(), isBox);
            }
        }
        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*temperature,    "temperature", isBox);
        writer.attachDofData(*phasePresence,  "phase presence", isBox);

        // if (velocityOutput.enableOutput()) // check if velocity output is demanded
        // {
        //     writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
        //     writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        // }

        writer.attachCellData(*rank, "process rank");
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one entity for the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void serializeEntity(std::ostream &outStream, const Entity &entity)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, entity);

        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize entity " << dofIdxGlobal);

        outStream << priVarSwitch().phasePresence(dofIdxGlobal) << " ";
    }

    /*!
     * \brief Reads the current solution from a restart file.
     *
     * \param inStream The input stream of one entity from the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, entity);

        // read phase presence
        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!inStream.good())
            DUNE_THROW(Dune::IOError, "Could not deserialize entity " << dofIdxGlobal);

        int phasePresence;
        inStream >> phasePresence;

        priVarSwitch_().setPhasePresence(dofIdxGlobal, phasePresence);
        priVarSwitch_().setOldPhasePresence(dofIdxGlobal, phasePresence);
    }

    const Dumux::TwoPTwoCPrimaryVariableSwitch<TypeTag>& priVarSwitch() const
    { return switch_; }

protected:

    Dumux::TwoPTwoCPrimaryVariableSwitch<TypeTag>& priVarSwitch_()
    { return switch_; }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        ParentType::applyInitialSolution_();

        // initialize the primary variable switch
        priVarSwitch_().init(this->problem_());
    }

    //! the class handling the primary variable switch
    TwoPTwoCPrimaryVariableSwitch<TypeTag> switch_;
    bool switchFlag_;
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
