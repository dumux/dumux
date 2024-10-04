#ifndef DUMUX_MULTIDOMAIN_TRANSIENT_MULTIBINARY_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_TRANSIENT_MULTIBINARY_COUPLINGMANAGER_HH

#include <memory>
#include <dune/common/hybridutilities.hh>
#include "multibinarycouplingmanager.hh"

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief Extension of the multi binary coupling manager that can be useful for the implementation
 *  of some coupling managers for transient problems using forward/backward Euler in time
 * \note The current implementation of the "previous solution" mechanism is a workaround that is
 * needed for handling the previous solution in the class `FCStaggeredFreeFlowCouplingManager`.
 * This approach is not optimal and should be refactored in future iterations.
 *
 * \tparam MDTraits the multidomain traits
 * \tparam CouplingMap a coupling policy class
 * \tparam CouplingMgrs the binary sub-coupling manager types
 */
template<class MDTraits, class CouplingMap, class ...CouplingMgrs>
class TransientMultiBinaryCouplingManager
: public MultiBinaryCouplingManager<MDTraits, CouplingMap, CouplingMgrs...>
{
    using ParentType = MultiBinaryCouplingManager<MDTraits, CouplingMap, CouplingMgrs...>;
    using SolutionVectors = typename ParentType::SolutionVectors;

public:
    TransientMultiBinaryCouplingManager()
    : MultiBinaryCouplingManager<MDTraits, CouplingMap, CouplingMgrs...>()
    {
        using namespace Dune::Hybrid;
        forEach(prevSolutionVectors_, [&](auto&& prevSolutionVector)
        {
            prevSolutionVector = std::make_shared<typename std::decay_t<decltype(prevSolutionVector)>::element_type>();
        });
    }

    //! Update the previous solution vector before assembly
    void updatePrevSolution(const typename MDTraits::SolutionVector& prevSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(prevSolutionVectors_)), [&](const auto id)
        {
            *std::get<id>(prevSolutionVectors_) = prevSol[id];
        });
    }

    SolutionVectors& prevSol()
    { return prevSolutionVectors_; }

    const SolutionVectors& prevSol() const
    { return prevSolutionVectors_; }

private:
    SolutionVectors prevSolutionVectors_;
};

} // end namespace Dumux

#endif
