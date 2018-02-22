// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup PorousmediumCompositional
 * \brief Reference implementation of a controller class for the Newton solver.
 *
 * Usually this controller should be sufficient.
 */
#ifndef DUMUX_PRIVARSWITCH_NEWTON_SOLVER_HH
#define DUMUX_PRIVARSWITCH_NEWTON_SOLVER_HH

#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux
{

/*!
 * \ingroup PorousmediumCompositional
 * \brief A newton controller that handles primary variable switches
 * \todo make this independent of TypeTag by making PrimaryVariableSwitch a template argument
 *       and extracting everything model specific from there
 * \todo Implement for volume variable caching enabled
 */
template <class TypeTag, class Assembler, class LinearSolver>
class PriVarSwitchNewtonSolver : public NewtonSolver<Assembler, LinearSolver>
{
    using Scalar =  typename Assembler::Scalar;
    using ParentType = NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename Assembler::ResidualType;
    using PrimaryVariableSwitch =  typename GET_PROP_TYPE(TypeTag, PrimaryVariableSwitch);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr auto discMethod = Assembler::FVGridGeometry::discMethod;
    static constexpr bool isBox = discMethod == DiscretizationMethod::box;

public:
    using ParentType::ParentType;

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const final
    {
        if (switchedInLastIteration_)
            return false;

        return ParentType::newtonConverged();
    }

    /*!
     *
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    void newtonBegin(const SolutionVector &u) final
    {
        ParentType::newtonBegin(u);
        priVarSwitch_ = std::make_unique<PrimaryVariableSwitch>(u.size());
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter) final
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        auto& assembler = this->assembler();
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        const auto& problem = assembler.problem();
        auto& gridVariables = assembler.gridVariables();

        // invoke the primary variable switch
        switchedInLastIteration_ = priVarSwitch_->update(uCurrentIter, gridVariables,
                                                         problem, fvGridGeometry);

        if(switchedInLastIteration_)
        {
            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                updateSwitchedVolVars_(std::integral_constant<bool, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>(),
                                       element, assembler, uCurrentIter, uLastIter);

                // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                // (not needed for box discretization)
                updateSwitchedFluxVarsCache_(std::integral_constant<bool, (GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache) && !isBox)>(),
                                             element, assembler, uCurrentIter, uLastIter);
            }
        }
    }

    /*!
     * \brief Called if the Newton method ended
     *        (not known yet if we failed or succeeded)
     */
    void newtonEnd() final
    {
        ParentType::newtonEnd();

        // in any way, we have to reset the switch flag
        switchedInLastIteration_ = false;
        // free some memory
        priVarSwitch_.release();
    }

private:

    /*!
     * \brief Update the volume variables whose primary variables were
              switched. Required when volume variables are cached globally.
     */
    template<class Element>
    void updateSwitchedVolVars_(std::true_type,
                                const Element& element,
                                Assembler& assembler,
                                const SolutionVector &uCurrentIter,
                                const SolutionVector &uLastIter)
    {
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        const auto& problem = assembler.problem();
        auto& gridVariables = assembler.gridVariables();

        // make sure FVElementGeometry is bound to the element
        auto fvGeometry = localView(fvGridGeometry);
        fvGeometry.bindElement(element);

        // update the secondary variables if global caching is enabled
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdxGlobal = scv.dofIndex();
            if (priVarSwitch_->wasSwitched(dofIdxGlobal))
            {
                const ElementSolutionVector elemSol(element, uCurrentIter, fvGridGeometry);
                auto& volVars = gridVariables.curGridVolVars().volVars(scv);
                volVars.update(elemSol, problem, element, scv);
            }
        }
    }

    /*!
     * \brief Update the fluxVars cache for dof whose primary variables were
              switched. Required when flux variables are cached globally.
     */
     template<class Element>
     void updateSwitchedFluxVarsCache_(std::true_type,
                                       const Element& element,
                                       Assembler& assembler,
                                       const SolutionVector &uCurrentIter,
                                       const SolutionVector &uLastIter)
    {
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        auto& gridVariables = assembler.gridVariables();

        // update the flux variables if global caching is enabled
        const auto dofIdxGlobal = fvGridGeometry.dofMapper().index(element);

        if (priVarSwitch_->wasSwitched(dofIdxGlobal))
        {
            // make sure FVElementGeometry and the volume variables are bound
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bind(element);
            auto curElemVolVars = localView(gridVariables.curGridVolVars());
            curElemVolVars.bind(element, fvGeometry, uCurrentIter);
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
        }
    }

     //! brief Do nothing when volume variables are not cached globally.
    template <typename... Args>
    void updateSwitchedVolVars_(std::false_type, Args&&... args) const {}

    //! brief Do nothing when flux variables are not cached globally.
    template <typename... Args>
    void updateSwitchedFluxVarsCache_(std::false_type, Args&&... args) const {}

    //! the class handling the primary variable switch
    std::unique_ptr<PrimaryVariableSwitch> priVarSwitch_;
    //! if we switched primary variables in the last iteration
    bool switchedInLastIteration_ = false;
};

} // end namespace Dumux

#endif
