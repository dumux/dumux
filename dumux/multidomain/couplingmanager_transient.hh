// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain, transient problems
 */

#ifndef DUMUX_MULTIDOMAIN_COUPLING_MANAGER_TRANSIENT_HH
#define DUMUX_MULTIDOMAIN_COUPLING_MANAGER_RANSIENT_HH

#include <memory>
#include <tuple>
#include <cstddef>
#include <utility>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/hybridutilities.hh>

#include "couplingmanager.hh"

namespace Dumux {


/*!
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for transient, multi domain problems,
 *  where staggered grid is coupled to another domain
 * \note The current implementation of the "previous solution" mechanism is a workaround that is
 * needed for the handling of the previous solution in the class `FCStaggeredFreeFlowCouplingManager`.
 * This approach is not optimal and should be refactored in future iterations.
 */
template<class Traits>
class TransientCouplingManager
: public CouplingManager<Traits>
{
    using ParentType = CouplingManager<Traits>;

protected:
    //! the type in which the solution vector is stored in the manager
    using SolutionVectorStorage = typename ParentType::SolutionVectorStorage;

public:
    /*!
     * \brief Default constructor
     *
     * The transient coupling manager stores the sub-solution vectors of the previous time-step, additional to the
     * sub-solution vectors stored in the parent coupling manager. For the non-owning pointer, attach the solution
     * vector managed elsewhere using `attachPrevSolution` and make sure that object stays alive of the lifetime
     * of the transient coupling manager.
     */
    TransientCouplingManager() : CouplingManager<Traits>()
    {
        using namespace Dune::Hybrid;
        forEach(prevSols_, [](auto& solutionVector){
            solutionVector = std::make_shared<typename std::decay_t<decltype(solutionVector)>::element_type>();
        });
    }

protected:
    /*!
     * \brief Attach a previous solution vector stored outside of this class.
     * \note The caller has to make sure that prevSol stays alive for the lifetime of
     *       the coupling manager. Otherwise we have a dangling reference here. Use with care.
     */
    void attachPrevSolution(typename ParentType::SolutionVectorStorage& prevSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(prevSols_)), [&](const auto id)
        {
            // do not take ownership of the external pointer's object
            std::get<id>(prevSols_) = Dune::stackobject_to_shared_ptr(*std::get<id>(prevSol));
        });
    }

    /*!
     * \brief the previous solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     * \return reference to the previous solution vector of the subproblem of type SubSolutionVector<i>
     */
    template<std::size_t i>
    auto& prevSol(Dune::index_constant<i> domainIdx)
    { return *std::get<i>(prevSols_); }

    /*!
     * \brief the previous solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     * \return const reference to the previous solution vector of the subproblem of type SubSolutionVector<i>
     */
    template<std::size_t i>
    const auto& prevSol(Dune::index_constant<i> domainIdx) const
    { return *std::get<i>(prevSols_); }

private:
    /*!
     * \brief A tuple of shared_ptr's to previous solution vectors of the subproblems
     */
    SolutionVectorStorage prevSols_;
};

} // end namespace Dumux

#endif
