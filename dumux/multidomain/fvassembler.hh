// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 *        with multiple domains
 */
#ifndef DUMUX_MULTIDOMAIN_FV_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_FV_ASSEMBLER_HH

#include <type_traits>
#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/common/gridcapabilities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/parallel/multithreading.hh>

#include "couplingjacobianpattern.hh"
#include "subdomaincclocalassembler.hh"
#include "subdomaincvfelocalassembler.hh"
#include "subdomainstaggeredlocalassembler.hh"
#include "subdomainfclocalassembler.hh"
#include "assemblerview.hh"

#include <dumux/discretization/method.hh>

namespace Dumux {

namespace Grid::Capabilities {

namespace Detail {
// helper for multi-domain models
template<class T, std::size_t... I>
bool allGridsSupportsMultithreadingImpl(const T& gridGeometries, std::index_sequence<I...>)
{
    return (... && supportsMultithreading(std::get<I>(gridGeometries)->gridView()));
}
} // end namespace Detail

// helper for multi-domain models (all grids have to support multithreading)
template<class... GG>
bool allGridsSupportsMultithreading(const std::tuple<GG...>& gridGeometries)
{
    return Detail::allGridsSupportsMultithreadingImpl<std::tuple<GG...>>(gridGeometries, std::make_index_sequence<sizeof...(GG)>());
}

} // end namespace Grid::Capabilities

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief trait that is specialized for coupling manager supporting multithreaded assembly
 */
template<class CM>
struct CouplingManagerSupportsMultithreadedAssembly : public std::false_type
{};

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 *        with multiple domains
 * \tparam MDTraits the multidimensional traits
 * \tparam diffMethod the differentiation method to residual compute derivatives
 * \tparam useImplicitAssembly if to use an implicit or explicit time discretization
 */
template<class MDTraits, class CMType, DiffMethod diffMethod, bool useImplicitAssembly = true>
class MultiDomainFVAssembler
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

public:
    using Traits = MDTraits;

    using Scalar = typename MDTraits::Scalar;

    //! TODO get rid of this GetPropType
    template<std::size_t id>
    using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id>
    using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;

    template<std::size_t id>
    using GridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;

    template<std::size_t id>
    using Problem = typename MDTraits::template SubDomain<id>::Problem;

    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = typename MDTraits::ResidualVector;

    using CouplingManager = CMType;

    /*!
     * \brief Returns true if the assembler considers implicit assembly.
     */
    static constexpr bool isImplicit()
    { return useImplicitAssembly; }

private:

    using ProblemTuple = typename MDTraits::template TupleOfSharedPtrConst<Problem>;
    using GridGeometryTuple = typename MDTraits::template TupleOfSharedPtrConst<GridGeometry>;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    using TimeLoop = TimeLoopBase<Scalar>;
    using ThisType = MultiDomainFVAssembler<MDTraits, CouplingManager, diffMethod, isImplicit()>;

    template<std::size_t id>
    using SubDomainAssemblerView = MultiDomainAssemblerSubDomainView<ThisType, id>;

    template<class DiscretizationMethod, std::size_t id>
    struct SubDomainAssemblerType;

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::CCTpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::CCMpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, isImplicit()>;
    };

    template<std::size_t id, class DM>
    struct SubDomainAssemblerType<DiscretizationMethods::CVFE<DM>, id>
    {
        using type = SubDomainCVFELocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::Staggered, id>
    {
        using type = SubDomainStaggeredLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::FCStaggered, id>
    {
        using type = SubDomainFaceCenteredLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    using SubDomainAssembler = typename SubDomainAssemblerType<typename GridGeometry<id>::DiscretizationMethod, id>::type;

public:


    /*!
     * \brief The constructor for stationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainFVAssembler(ProblemTuple problem,
                           GridGeometryTuple gridGeometry,
                           GridVariablesTuple gridVariables,
                           std::shared_ptr<CouplingManager> couplingManager)
    : couplingManager_(couplingManager)
    , problemTuple_(std::move(problem))
    , gridGeometryTuple_(std::move(gridGeometry))
    , gridVariablesTuple_(std::move(gridVariables))
    , timeLoop_()
    , isStationaryProblem_(true)
    , warningIssued_(false)
    {
        static_assert(isImplicit(), "Explicit assembler for stationary problem doesn't make sense!");
        std::cout << "Instantiated assembler for a stationary problem." << std::endl;

        enableMultithreading_ = CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value
            && Grid::Capabilities::allGridsSupportsMultithreading(gridGeometryTuple_)
            && !Multithreading::isSerial()
            && getParam<bool>("Assembly.Multithreading", true);

        maybeComputeColors_();
    }

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainFVAssembler(ProblemTuple problem,
                           GridGeometryTuple gridGeometry,
                           GridVariablesTuple gridVariables,
                           std::shared_ptr<CouplingManager> couplingManager,
                           std::shared_ptr<const TimeLoop> timeLoop,
                           const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , problemTuple_(std::move(problem))
    , gridGeometryTuple_(std::move(gridGeometry))
    , gridVariablesTuple_(std::move(gridVariables))
    , timeLoop_(timeLoop)
    , prevSol_(&prevSol)
    , isStationaryProblem_(false)
    , warningIssued_(false)
    {
        std::cout << "Instantiated assembler for an instationary problem." << std::endl;

        enableMultithreading_ = CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value
            && Grid::Capabilities::allGridsSupportsMultithreading(gridGeometryTuple_)
            && !Multithreading::isSerial()
            && getParam<bool>("Assembly.Multithreading", true);

        maybeComputeColors_();
    }

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        checkAssemblerState_();
        resetJacobian_();
        resetResidual_();

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& jacRow = (*jacobian_)[domainId];
            auto& subRes = (*residual_)[domainId];
            this->assembleJacobianAndResidual_(domainId, jacRow, subRes, curSol);

            const auto gridGeometry = std::get<domainId>(gridGeometryTuple_);
            enforcePeriodicConstraints_(domainId, jacRow, subRes, *gridGeometry, curSol[domainId]);
        });
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        assembleResidual(*residual_, curSol);
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const SolutionVector& curSol)
    {
        r = 0.0;

        checkAssemblerState_();

        // update the grid variables for the case of active caching
        updateGridVariables(curSol);

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(r)), [&](const auto domainId)
        {
            auto& subRes = r[domainId];
            this->assembleResidual_(domainId, subRes, curSol);
        });
    }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<ResidualType> r)
    {
        jacobian_ = A;
        residual_ = r;

        setJacobianBuildMode(*jacobian_);
        setJacobianPattern_(*jacobian_);
        setResidualSize_(*residual_);
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        residual_ = std::make_shared<ResidualType>();

        setJacobianBuildMode(*jacobian_);
        setJacobianPattern_(*jacobian_);
        setResidualSize_(*residual_);
    }

    /*!
     * \brief Sets the jacobian build mode
     */
    void setJacobianBuildMode(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto i)
        {
            forEach(jac[i], [&](auto& jacBlock)
            {
                using BlockType = std::decay_t<decltype(jacBlock)>;
                if (jacBlock.buildMode() == BlockType::BuildMode::unknown)
                    jacBlock.setBuildMode(BlockType::BuildMode::random);
                else if (jacBlock.buildMode() != BlockType::BuildMode::random)
                    DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");
            });
        });
    }

     /*!
     * \brief Resizes jacobian and residual and recomputes colors
     */
    void updateAfterGridAdaption()
    {
        setJacobianPattern_(*jacobian_);
        setResidualSize_(*residual_);
        maybeComputeColors_();
    }

    /*!
     * \brief Updates the grid variables with the given solution
     */
    void updateGridVariables(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).update(curSol[domainId]); });
    }

    /*!
     * \brief Resets the grid variables to the last time step
     */
    void resetTimeStep(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).resetTimeStep(curSol[domainId]); });
    }

    //! the number of dof locations of domain i
    template<std::size_t i>
    std::size_t numDofs(Dune::index_constant<i> domainId) const
    { return std::get<domainId>(gridGeometryTuple_)->numDofs(); }

    //! the problem of domain i
    template<std::size_t i>
    const auto& problem(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(problemTuple_); }

    //! the finite volume grid geometry of domain i
    template<std::size_t i>
    const auto& gridGeometry(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridGeometryTuple_); }

    //! the grid view of domain i
    template<std::size_t i>
    const auto& gridView(Dune::index_constant<i> domainId) const
    { return gridGeometry(domainId).gridView(); }

    //! the grid variables of domain i
    template<std::size_t i>
    GridVariables<i>& gridVariables(Dune::index_constant<i> domainId)
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the grid variables of domain i
    template<std::size_t i>
    const GridVariables<i>& gridVariables(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    //! the full Jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! the full residual vector
    ResidualType& residual()
    { return *residual_; }

    //! the solution of the previous time step
    const SolutionVector& prevSol() const
    { return *prevSol_; }

    /*!
     * \brief Set time loop for instationary problems
     * \note calling this turns this into a stationary assembler
     */
    void setTimeManager(std::shared_ptr<const TimeLoop> timeLoop)
    { timeLoop_ = timeLoop; isStationaryProblem_ = !(static_cast<bool>(timeLoop)); }

    /*!
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { prevSol_ = &u; }

    /*!
     * \brief Whether we are assembling a stationary or instationary problem
     */
    bool isStationaryProblem() const
    { return isStationaryProblem_; }

    /*!
     * \brief Create a local residual object (used by the local assembler)
     */
    template<std::size_t i>
    LocalResidual<i> localResidual(Dune::index_constant<i> domainId) const
    { return LocalResidual<i>(std::get<domainId>(problemTuple_).get(), timeLoop_.get()); }

protected:
    //! the coupling manager coupling the sub domains
    std::shared_ptr<CouplingManager> couplingManager_;

private:
    /*!
     * \brief Sets the jacobian sparsity pattern.
     */
    void setJacobianPattern_(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(jac[domainI])), [&](const auto domainJ)
            {
                const auto pattern = this->getJacobianPattern_(domainI, domainJ);
                pattern.exportIdx(jac[domainI][domainJ]);
            });
        });
    }

    /*!
     * \brief Resizes the residual
     */
    void setResidualSize_(ResidualType& res) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(res)), [&](const auto domainId)
        { res[domainId].resize(this->numDofs(domainId)); });
    }

    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if(!residual_)
        {
            residual_ = std::make_shared<ResidualType>();
            setResidualSize_(*residual_);
        }

        (*residual_) = 0.0;
    }

    // reset the jacobian vector to 0.0
    void resetJacobian_()
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            setJacobianBuildMode(*jacobian_);
            setJacobianPattern_(*jacobian_);
        }

       (*jacobian_)  = 0.0;
    }

    //! Computes the colors
    void maybeComputeColors_()
    {
        if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            if (enableMultithreading_)
                couplingManager_->computeColorsForAssembly();
    }

    // check if the assembler is in a correct state for assembly
    void checkAssemblerState_() const
    {
        if (!isStationaryProblem_ && !prevSol_)
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        if (isStationaryProblem_ && prevSol_)
            DUNE_THROW(Dune::InvalidStateException, "Assembling stationary problem but a previous solution was set."
                                                    << " Did you forget to set the timeLoop to make this problem instationary?");
    }

    template<std::size_t i, class JacRow, class SubRes>
    void assembleJacobianAndResidual_(Dune::index_constant<i> domainId, JacRow& jacRow, SubRes& subRes,
                                      const SolutionVector& curSol)
    {
        assemble_(domainId, [&](const auto& element)
        {
            MultiDomainAssemblerSubDomainView view{*this, domainId};
            SubDomainAssembler<i> subDomainAssembler(view, element, curSol, *couplingManager_);
            subDomainAssembler.assembleJacobianAndResidual(jacRow, subRes, gridVariablesTuple_);
        });
    }

    template<std::size_t i, class SubRes>
    void assembleResidual_(Dune::index_constant<i> domainId, SubRes& subRes,
                           const SolutionVector& curSol)
    {
        assemble_(domainId, [&](const auto& element)
        {
            MultiDomainAssemblerSubDomainView view{*this, domainId};
            SubDomainAssembler<i> subDomainAssembler(view, element, curSol, *couplingManager_);
            subDomainAssembler.assembleResidual(subRes);
        });
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if an exception is thrown during assembly
     */
    template<std::size_t i, class AssembleElementFunc>
    void assemble_(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
        try
        {
            if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            {
                if (enableMultithreading_)
                {
                    couplingManager_->assembleMultithreaded(
                        domainId, std::forward<AssembleElementFunc>(assembleElement)
                    );
                    return;
                }
            }

            // fallback for coupling managers that don't support multithreaded assembly (yet)
            // or if multithreaded assembly is disabled
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView(domainId)))
                assembleElement(element);

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem occurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView(domainId).comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView(domainId).comm().size() > 1)
            succeeded = gridView(domainId).comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    // get diagonal block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i==j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        const auto& gg = gridGeometry(domainI);
        auto pattern = getJacobianPattern<isImplicit()>(gg);
        couplingManager_->extendJacobianPattern(domainI, pattern);
        return pattern;
    }

    // get coupling block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i!=j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        return getCouplingJacobianPattern<isImplicit()>(*couplingManager_,
                                                        domainI, gridGeometry(domainI),
                                                        domainJ, gridGeometry(domainJ));
    }

    // build periodic constraints into the system matrix
    template<std::size_t i, class JacRow, class Res, class GG, class Sol>
    void enforcePeriodicConstraints_(Dune::index_constant<i> domainI, JacRow& jacRow, Res& res, const GG& gridGeometry, const Sol& curSol)
    {
        if constexpr (GG::discMethod == DiscretizationMethods::box || GG::discMethod == DiscretizationMethods::fcstaggered)
        {
            for (const auto& m : gridGeometry.periodicVertexMap())
            {
                if (m.first < m.second)
                {
                    auto& jac = jacRow[domainI];

                    // add the second row to the first
                    res[m.first] += res[m.second];

                    // enforce the solution of the first periodic DOF to the second one
                    res[m.second] = curSol[m.second] - curSol[m.first];

                    const auto end = jac[m.second].end();
                    for (auto it = jac[m.second].begin(); it != end; ++it)
                        jac[m.first][it.index()] += (*it);

                    // enforce constraint in second row
                    for (auto it = jac[m.second].begin(); it != end; ++it)
                        (*it) = it.index() == m.second ? 1.0 : it.index() == m.first ? -1.0 : 0.0;

                    using namespace Dune::Hybrid;
                    forEach(makeIncompleteIntegerSequence<JacRow::size(), domainI>(), [&](const auto couplingDomainId)
                    {
                        auto& jacCoupling = jacRow[couplingDomainId];

                        for (auto it = jacCoupling[m.second].begin(); it != jacCoupling[m.second].end(); ++it)
                            jacCoupling[m.first][it.index()] += (*it);

                        for (auto it = jacCoupling[m.second].begin(); it != jacCoupling[m.second].end(); ++it)
                            (*it) = 0.0;
                    });
                }
            }
        }
    }

    //! pointer to the problem to be solved
    ProblemTuple problemTuple_;

    //! the finite volume geometry of the grid
    GridGeometryTuple gridGeometryTuple_;

    //! the variables container for the grid
    GridVariablesTuple gridVariablesTuple_;

    //! the time loop for instationary problem assembly
    std::shared_ptr<const TimeLoop> timeLoop_;

    //! an observing pointer to the previous solution for instationary problems
    const SolutionVector* prevSol_ = nullptr;

    //! if this assembler is assembling an instationary problem
    bool isStationaryProblem_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualType> residual_;

    //! Issue a warning if the calculation is used in parallel with overlap. This could be a static local variable if it wasn't for g++7 yielding a linker error.
    mutable bool warningIssued_;

    //! if multithreaded assembly is enabled
    bool enableMultithreading_ = false;
};

} // end namespace Dumux

#endif
