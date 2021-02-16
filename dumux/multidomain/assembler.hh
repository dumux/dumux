// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 *        with multiple domains
 */
#ifndef DUMUX_MULTIDOMAIN_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_ASSEMBLER_HH

#include <type_traits>
#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/jacobianpattern.hh>

#include <dumux/multidomain/fvgridvariables.hh>
#include <dumux/multidomain/fvgridgeometry.hh>
#include <dumux/timestepping/multistagetimestepper.hh>

#include "couplingjacobianpattern.hh"
// #include "subdomaincclocalassembler.hh"
// #include "subdomainboxlocalassembler.hh"
#include "subdomainstaggeredlocalassembler.hh"

#include <dumux/discretization/method.hh>
#include "assembly/localassembler.hh"

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 *        with multiple domains
 * \tparam MDTraits the multidimension traits
 * \tparam diffMethod the differentiation method to residual compute derivatives
 * \tparam useImplicitAssembly if to use an implicit or explicit time discretization
 */
template<class MDTraits, class CMType, DiffMethod diffMethod, bool useImplicitAssembly = true>
class MultiDomainAssembler
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

public:
    using Traits = MDTraits;
    using CouplingManager = CMType;

    // types related to linear algebra
    using Scalar = typename MDTraits::Scalar;
    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = SolutionVector;

    // types underlying the subdomains
    template<std::size_t id> using SubDomainLocalOperator = typename MDTraits::template SubDomain<id>::LocalOperator;
    template<std::size_t id> using SubDomainGridVariables = typename MDTraits::template SubDomain<id>::GridVariables;
    template<std::size_t id> using SubDomainGridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;
    template<std::size_t id> using SubDomainProblem = typename MDTraits::template SubDomain<id>::Problem;

    // The variables that define an evaluation point
    using Variables = MultiDomainFVGridVariables<Traits>;
    using GridVariables = Variables;

    // The grid geometry on which it is operated
    using GridGeometry = MultiDomainFVGridGeometry<Traits>;

    //! export the type for parameters of a stage in time integration
    using StageParams = MultiStageParams<Scalar>;

private:
    template<std::size_t id> using SubDomainGridView = typename SubDomainGridGeometry<id>::GridView;
    template<std::size_t id> using SubDomainElemVars = typename SubDomainGridVariables<id>::LocalView;
    template<std::size_t id> using SubDomainElement = typename SubDomainGridView<id>::template Codim<0>::Entity;

    using TimeLoop = TimeLoopBase<Scalar>;
    using ThisType = MultiDomainAssembler<MDTraits, CouplingManager, diffMethod>;
    using GridGeometryTuple = typename MDTraits::template TupleOfSharedPtrConst<SubDomainGridGeometry>;

    template<std::size_t id>
    using Context = typename CouplingManager::template CouplingContext<id>;

public:

    /*!
     * \brief Constructor from a grid geometry tuple for stationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainAssembler(GridGeometryTuple&& gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : gridGeometryTuple_(gridGeometry)
    , couplingManager_(couplingManager)
    , isImplicit_(true)
    {
        // TODO: Store that instantiation happened for stationary problem for later error catch?
        //       Or maybe we need another mechanism to tell the assembler if implicit/explicit...
        std::cout << "Instantiated MultiDomainAssembler for a stationary problem." << std::endl;
    }

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    template<class TimeIntegrationMethod>
    MultiDomainAssembler(GridGeometryTuple&& gridGeometry,
                         std::shared_ptr<CouplingManager> couplingManager,
                         const TimeIntegrationMethod& method)
    : gridGeometryTuple_(gridGeometry)
    , couplingManager_(couplingManager)
    , isImplicit_(method.implicit())
    {
        std::cout << "Instantiated MultiDomainAssembler for an instationary problem." << std::endl;
    }

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const GridVariables& gridVars)
    {
        ptr_ = &gridVars;
        resetJacobian_();
        resetResidual_();

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& jacRow = (*jacobian_)[domainId];
            auto& subRes = (*residual_)[domainId];
            this->assembleJacobianAndResidual_(domainId, jacRow, subRes, gridVars);
        });

        using namespace Dune::Indices;
        // std::cout << "SOL in GV: " << std::endl;
        // Dune::printvector(std::cout, gridVars.dofs()[_0], "", "", 1, 10, 5);
        // Dune::printvector(std::cout, gridVars.dofs()[_1], "", "", 1, 10, 5);
        //
        // std::cout << "Linear system: " << std::endl;
        // Dune::printmatrix(std::cout, (*jacobian_)[_0][_0], "", "", 10, 2);
        // Dune::printmatrix(std::cout, (*jacobian_)[_0][_1], "", "", 10, 2);
        // Dune::printmatrix(std::cout, (*jacobian_)[_1][_0], "", "", 10, 2);
        // Dune::printmatrix(std::cout, (*jacobian_)[_1][_1], "", "", 10, 2);
        // Dune::printvector(std::cout, (*residual_)[_0], "", "", 1, 10, 2);
        // Dune::printvector(std::cout, (*residual_)[_1], "", "", 1, 10, 2);
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const GridVariables& gridVars)
    {
        ptr_ = &gridVars;
        resetResidual_();
        assembleResidual(*residual_, gridVars);
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const GridVariables& gridVars)
    {
        ptr_ = &gridVars;
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(r)), [&](const auto domainId)
        {
            auto& subRes = r[domainId];
            this->assembleResidual_(domainId, subRes, gridVars);
        });
    }

    //! compute the residual and return it's vector norm
    //! TODO: this needs to be adapted in parallel
    Scalar residualNorm(const GridVariables& gridVars)
    {
        ptr_ = &gridVars;
        ResidualType residual;
        setResidualSize(residual);
        assembleResidual(residual, gridVars);

        // calculate the square norm of the residual
        const Scalar result2 = residual.two_norm2();
        using std::sqrt;
        return sqrt(result2);
    }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        setJacobianBuildMode(*jacobian_);
        setJacobianPattern(*jacobian_);
        setResidualSize(*residual_);
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        residual_ = std::make_shared<SolutionVector>();

        setJacobianBuildMode(*jacobian_);
        setJacobianPattern(*jacobian_);
        setResidualSize(*residual_);
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
     * \brief Sets the jacobian sparsity pattern.
     */
    void setJacobianPattern(JacobianMatrix& jac) const
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
    void setResidualSize(SolutionVector& res) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(res)), [&](const auto domainId)
        { res[domainId].resize(this->numDofs(domainId)); });
    }

    //! the number of dof locations of domain i
    template<std::size_t i>
    std::size_t numDofs(Dune::index_constant<i> domainId) const
    { return std::get<domainId>(gridGeometryTuple_)->numDofs(); }

    //! the problem of domain i
    // template<std::size_t i>
    // const auto& problem(Dune::index_constant<i> domainId) const
    // { return *std::get<domainId>(problemTuple_); }

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
    auto& gridVariables(Dune::index_constant<i> domainId)
    { return (*ptr_)[domainId]; }

    //! the grid variables of domain i
    template<std::size_t i>
    const auto& gridVariables(Dune::index_constant<i> domainId) const
    { return (*ptr_)[domainId]; }

    //! the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    //! the full Jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! the full residual vector
    SolutionVector& residual()
    { return *residual_; }

    /*!
     * \brief Prepare for a new stage within a time integration step.
     *        This caches the given grid variables, which are then used as a
     *        representation of the previous stage. Moreover, the given grid
     *        variables are then updated to the time level of the upcoming stage.
     * \param gridVars the grid variables representing the previous stage
     * \param params the parameters with the weights to be used in the upcoming stage
     * \todo TODO: This function does two things, namely caching and then updating.
     *             Should we split/delegate this, or is the current name descriptive enough?
     *             When used from outside, one would expect the gridvars to be prepared maybe,
     *             and that is what's done. Caching might not be expected from the outside but
     *             it is also not important that that is known from there?
     */
    void prepareStage(GridVariables& gridVars,
                      std::shared_ptr<const StageParams> params)
    {
        stageParams_ = params;
        const auto curStage = params->size() - 1;

        // we keep track of previous stages, they are needed for residual assembly
        prevStageVariables_.push_back(gridVars.deepCopy());

        // Now we update the time level of the given grid variables
        const auto t = params->timeAtStage(curStage);
        const auto prevT = params->timeAtStage(0);
        const auto dtFraction = params->timeStepFraction(curStage);
        TimeLevel<Scalar> timeLevel(t, prevT, dtFraction);

        gridVars.updateTime(timeLevel);
    }

    /*!
     * \brief Remove traces from stages within a time integration step.
     */
    void clearStages()
    {
        prevStageVariables_.clear();
        stageParams_ = nullptr;
    }

private:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize(*residual_);
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
            setJacobianPattern(*jacobian_);
        }

       (*jacobian_)  = 0.0;
    }

    template<std::size_t i, class JacRow, class SubRes>
    void assembleJacobianAndResidual_(Dune::index_constant<i> domainId,
                                      JacRow& jacRow, SubRes& subRes,
                                      const GridVariables& gridVars)
    {
        assemble_(domainId, [&](const auto& element)
        {
            auto ggLocalView = localView(gridGeometry(domainId));
            ggLocalView.bind(element);
            auto [elemVars, contexts] = this->prepareElemVariables_<domainId>(gridVars, element, ggLocalView);

            using LocalAssembler = Dumux::MultiDomainLocalAssembler<i, ThisType, diffMethod>;
            LocalAssembler localAssembler(element, ggLocalView, contexts, elemVars, stageParams_, couplingManager_);
            localAssembler.assembleJacobianAndResidual(jacRow, subRes);
        });
    }

    template<std::size_t i, class SubRes>
    void assembleResidual_(Dune::index_constant<i> domainId,
                           SubRes& subRes,
                           const GridVariables& gridVars)
    {
        assemble_(domainId, [&](const auto& element)
        {
            auto ggLocalView = localView(gridGeometry(domainId));
            ggLocalView.bind(element);
            auto [elemVars, contexts] = this->prepareElemVariables_<domainId>(gridVars, element, ggLocalView);

            using LocalAssembler = Dumux::MultiDomainLocalAssembler<i, ThisType, diffMethod>;
            LocalAssembler localAssembler(element, ggLocalView, contexts, elemVars, stageParams_, couplingManager_);
            localAssembler.assembleResidual(subRes);
        });
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     * TODO: assemble in parallel
     */
    template<std::size_t i, class AssembleElementFunc>
    void assemble_(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        // let the local assembler add the element contributions
        for (const auto& element : elements(gridView(domainId)))
            assembleElement(element);
    }

    // get diagonal block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i==j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        const auto& gg = gridGeometry(domainI);
        auto pattern = isImplicit_ ? getJacobianPattern<true>(gg)
                                   : getJacobianPattern<false>(gg);
        couplingManager_->extendJacobianPattern(domainI, pattern);
        return pattern;
    }

    // get coupling block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i!=j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        return isImplicit_ ? getCouplingJacobianPattern<true>(*couplingManager_,
                                                              domainI, gridGeometry(domainI),
                                                              domainJ, gridGeometry(domainJ))
                           : getCouplingJacobianPattern<false>(*couplingManager_,
                                                               domainI, gridGeometry(domainI),
                                                               domainJ, gridGeometry(domainJ));
    }

    //! prepares the local views on the grid variables for the given element
    //! \todo: TODO: when stageparams.skipSpatial() == true, we don't need to bind flux vars caches!
    template<std::size_t id>
    std::pair< std::vector<SubDomainElemVars<id>>,
               std::vector<std::shared_ptr<Context<id>>> >
    prepareElemVariables_(const Variables& gridVars,
                          const SubDomainElement<id>& element,
                          const typename SubDomainGridGeometry<id>::LocalView& ggLocalView) const
    {
        static constexpr auto domainId = Dune::index_constant<id>();

        if (!stageParams_)
        {
            auto context = couplingManager_->makeCouplingContext(domainId, element, gridVars);
            auto elemVars = localView(gridVars[domainId]);
            elemVars.bind(element, ggLocalView);
            return std::make_pair(std::vector<SubDomainElemVars<id>>{{elemVars}},
                                  std::vector<std::shared_ptr<Context<id>>>{{context}});
        }
        else
        {
            std::vector<SubDomainElemVars<id>> elemVars;
            std::vector<std::shared_ptr<Context<id>>> contexts;

            elemVars.reserve(stageParams_->size());
            contexts.reserve(stageParams_->size());

            for (int i = 0; i < stageParams_->size()-1; ++i)
            {
                elemVars.emplace_back(prevStageVariables_[i][domainId]);
                contexts.emplace_back(couplingManager_->makeCouplingContext(domainId, element, prevStageVariables_[i]));
            }

            elemVars.emplace_back(gridVars[domainId]);
            contexts.emplace_back(couplingManager_->makeCouplingContext(domainId, element, gridVars));

            for (int i = 0; i < elemVars.size(); ++i)
            {
                couplingManager_->setCouplingContext(contexts[i]);
                elemVars[i].bind(element, ggLocalView);
            }

            return std::make_pair(std::vector<SubDomainElemVars<id>>{elemVars},
                                  std::vector<std::shared_ptr<Context<id>>>{contexts});
        }
    }

    //! the grid geometries of all subdomains
    GridGeometryTuple gridGeometryTuple_;

    //! pointer to the coupling manager
    std::shared_ptr<CouplingManager> couplingManager_;

    //! states if an implicit of explicit scheme is used (affects jacobian pattern)
    bool isImplicit_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<SolutionVector> residual_;

    //! parameters containing information on the current stage of time integration
    std::shared_ptr<const StageParams> stageParams_ = nullptr;

    //! keep track of the states of previous stages within a time integration step
    std::vector<GridVariables> prevStageVariables_;

    const GridVariables* ptr_;
};

} // end namespace Dumux

#endif
