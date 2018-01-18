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
#ifndef DUMUX_PRIVARSWITCH_NEWTON_CONTROLLER_HH
#define DUMUX_PRIVARSWITCH_NEWTON_CONTROLLER_HH

#include <memory>
#include <dune/common/hybridutilities.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux
{

/*!
 * \ingroup PorousmediumCompositional
 * \brief A newton controller that handles primary variable switches
 * \todo make this independent of TypeTag by making PrimaryVariableSwitch a template argument
 *       and extracting everything model specific from there
 * \todo Implement for volume variable caching enabled
 */
template <class TypeTag>
class PriVarSwitchNewtonController : public NewtonController<typename GET_PROP_TYPE(TypeTag, Scalar)>
{
    using Scalar =  typename GET_PROP_TYPE(TypeTag, Scalar);
    using ParentType = NewtonController<Scalar>;
    using PrimaryVariableSwitch =  typename GET_PROP_TYPE(TypeTag, PrimaryVariableSwitch);

    // using ElementSolution =  typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    // static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box;

public:
    using ParentType::ParentType;

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
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
    template<class SolutionVector>
    void newtonBegin(const SolutionVector &u)
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
    template<class JacobianAssembler, class SolutionVector>
    void newtonEndStep(JacobianAssembler& assembler,
                       SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ParentType::newtonEndStep(assembler, uCurrentIter, uLastIter);

        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        const auto& problem = assembler.problem();
        auto& gridVariables = assembler.gridVariables();

        switchedInLastIteration_ = priVarSwitch_->update(uCurrentIter, gridVariables,
                                                         problem, fvGridGeometry);

        Dune::Hybrid::ifElse(std::integral_constant<bool, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>(),
        [&](auto&& _if)
        {
            DUNE_THROW(Dune::NotImplemented, "Privar switch and volume varaible caching! Please implement!");

            // TODO:
            // std::cout << "blub";
            //
            // // update the secondary variables if global caching is enabled
            // // \note we only updated if phase presence changed as the volume variables
            // //       are already updated once by the switch
            // for (const auto& element : elements(fvGridGeometry.gridView()))
            // {
            //     // make sure FVElementGeometry & vol vars are bound to the element
            //     auto fvGeometry = localView(fvGridGeometry);
            //     fvGeometry.bindElement(element);
            //
            //     if (switchedInLastIteration_)
            //     {
            //         for (auto&& scv : scvs(fvGeometry))
            //         {
            //             const auto dofIdxGlobal = scv.dofIndex();
            //             if (priVarSwitch_->wasSwitched(dofIdxGlobal))
            //             {
            //                 const auto eIdx = fvGridGeometry.elementMapper().index(element);
            //                 const ElementSolution elemSol(element, this->curSol(), fvGridGeometry);
            //                 this->nonConstCurGridVolVars().volVars(eIdx, scv.indexInElement()).update(elemSol,
            //                                                                                             problem,
            //                                                                                             element,
            //                                                                                             scv);
            //             }
            //         }
            //     }
            //
            //     // handle the boundary volume variables for cell-centered models
            //     if(!isBox)
            //     {
            //         for (auto&& scvf : scvfs(fvGeometry))
            //         {
            //             // if we are not on a boundary, skip the rest
            //             if (!scvf.boundary())
            //                 continue;
            //
            //             // check if boundary is a pure dirichlet boundary
            //             const auto bcTypes = problem.boundaryTypes(element, scvf);
            //             if (bcTypes.hasOnlyDirichlet())
            //             {
            //                 const auto insideScvIdx = scvf.insideScvIdx();
            //                 const auto& insideScv = fvGeometry.scv(insideScvIdx);
            //                 const ElementSolution elemSol(problem.dirichlet(element, scvf));
            //
            //                 this->nonConstCurGridVolVars().volVars(scvf.outsideScvIdx(), 0/*indexInElement*/).update(elemSol, problem, element, insideScv);
            //             }
            //         }
            //     }
            // }
        });
    }

    /*!
     * \brief Called if the Newton method ended
     *        (not known yet if we failed or succeeded)
     */
    void newtonEnd()
    {
        ParentType::newtonEnd();

        // in any way, we have to reset the switch flag
        switchedInLastIteration_ = false;
        // free some memory
        priVarSwitch_.release();
    }

private:

    //! the class handling the primary variable switch
    std::unique_ptr<PrimaryVariableSwitch> priVarSwitch_;
    //! if we switched primary variables in the last iteration
    bool switchedInLastIteration_ = false;
};

} // end namespace Dumux

#endif
