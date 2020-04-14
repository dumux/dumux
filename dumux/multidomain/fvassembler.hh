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
#ifndef DUMUX_MULTIDOMAIN_FV_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_FV_ASSEMBLER_HH

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

#include "couplingjacobianpattern.hh"
#include "subdomaincclocalassembler.hh"
#include "subdomainboxlocalassembler.hh"
#include "subdomainstaggeredlocalassembler.hh"

#include <dumux/discretization/method.hh>

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
    using IndexType = typename GridGeometry<id>::GridView::IndexSet::IndexType;

    template<std::size_t id>
    using Problem = typename MDTraits::template SubDomain<id>::Problem;

    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = SolutionVector;

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

    template<DiscretizationMethod discMethod, std::size_t id>
    struct SubDomainAssemblerType;

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::cctpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::ccmpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::box, id>
    {
        using type = SubDomainBoxLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::staggered, id>
    {
        using type = SubDomainStaggeredLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    using SubDomainAssembler = typename SubDomainAssemblerType<GridGeometry<id>::discMethod, id>::type;

public:


    /*!
     * \brief The constructor for stationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainFVAssembler(ProblemTuple&& problem,
                           GridGeometryTuple&& gridGeometry,
                           GridVariablesTuple&& gridVariables,
                           std::shared_ptr<CouplingManager> couplingManager)
    : couplingManager_(couplingManager)
    , problemTuple_(problem)
    , gridGeometryTuple_(gridGeometry)
    , gridVariablesTuple_(gridVariables)
    , timeLoop_()
    , isStationaryProblem_(true)
    {
        static_assert(isImplicit(), "Explicit assembler for stationary problem doesn't make sense!");
        std::cout << "Instantiated assembler for a stationary problem." << std::endl;
    }

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainFVAssembler(ProblemTuple&& problem,
                           GridGeometryTuple&& gridGeometry,
                           GridVariablesTuple&& gridVariables,
                           std::shared_ptr<CouplingManager> couplingManager,
                           std::shared_ptr<const TimeLoop> timeLoop,
                           const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , problemTuple_(problem)
    , gridGeometryTuple_(gridGeometry)
    , gridVariablesTuple_(gridVariables)
    , timeLoop_(timeLoop)
    , prevSol_(&prevSol)
    , isStationaryProblem_(false)
    {
        std::cout << "Instantiated assembler for an instationary problem." << std::endl;
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
        forEach(integralRange(Dune::Hybrid::size(*jacobian_)), [&](const auto domainId)
        {
            auto& jacRow = (*jacobian_)[domainId];
            auto& subRes = (*residual_)[domainId];
            this->assembleJacobianAndResidual_(domainId, jacRow, subRes, curSol);
        });

        resetReducedCoefficientMatrix_();
        resetReducedRHS_();

        fillReducedCoefficientMatrix_(*jacobian_, *reducedCoefficientMatrix_, *furtherReducedCoefficientMatrix_);
        fillReducedRHS_(*residual_, *reducedResidual_);
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        assembleResidual(*residual_, curSol);

        resetReducedResidual_();
        fillReducedResidual_(*residual_, *reducedResidual_);
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const SolutionVector& curSol)
    {
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

    //! compute the residual and return it's vector norm
    //! TODO: this needs to be adapted in parallel
    Scalar residualNorm(const SolutionVector& curSol)
    {
        ResidualType residual;
        setResidualSize(residual);
        assembleResidual(residual, curSol);

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
     *        the coefficient matrix and right-hand side objects in this assembler if you don't need them outside this class
     */
    void setReducedLinearSystem()
    {
        reducedCoefficientMatrix_ = std::make_shared<JacobianMatrix>();
        furtherReducedCoefficientMatrix_ = std::make_shared<JacobianMatrix>();
        reducedResidual_ = std::make_shared<SolutionVector>();

        setJacobianBuildMode(*reducedCoefficientMatrix_);
        setJacobianBuildMode(*furtherReducedCoefficientMatrix_);

        setReducedCoefficientMatrixPattern(*jacobian_, *reducedCoefficientMatrix_, *furtherReducedCoefficientMatrix_);
        setReducedRHSSize(*reducedResidual_);
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
        forEach(jac, [](auto& jacRow)
        {
            forEach(jacRow, [](auto& jacBlock)
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
        forEach(integralRange(Dune::Hybrid::size(jac)), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(jac[domainI])), [&](const auto domainJ)
            {
                const auto pattern = this->getJacobianPattern_(domainI, domainJ);
                pattern.exportIdx(jac[domainI][domainJ]);
            });
        });
    }

    /*!
     * \brief Sets the reduced coefficient matrix' sparsity pattern.
     */
    void setReducedCoefficientMatrixPattern(const JacobianMatrix& coefficientMatrix, JacobianMatrix& reducedCoefficientMatrix, JacobianMatrix& furtherReducedCoefficientMatrix)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(coefficientMatrix)), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(coefficientMatrix[domainI])), [&](const auto domainJ)
            {
                // determine the dofs that do not take part in intersections
                auto& matrix = coefficientMatrix[domainI][domainJ];
                auto& reducedMatrix = reducedCoefficientMatrix[domainI][domainJ];

                std::vector<bool> isVoidRow(matrix.N(), false);
                const auto& v=(this->reductionIndexSet(domainI));

                for (std::size_t dofIdxTarget = 0; dofIdxTarget < matrix.N(); ++dofIdxTarget)
                {
                    if (std::find(v.begin(), v.end(), dofIdxTarget) != v.end())
                        isVoidRow[dofIdxTarget] = true;
                }

                //TODO deal with empty isVoidRow or isVoidCol
                this->setPatternDeleteSetOfRowsFromBCRSMatrix_(matrix, reducedMatrix, isVoidRow);

                const auto& w=(this->reductionIndexSet(domainJ));
                std::vector<bool> isVoidCol(reducedMatrix.M(), false);
                for (std::size_t dofIdxTarget = 0; dofIdxTarget < reducedMatrix.M(); ++dofIdxTarget)
                {
                    if (std::find(w.begin(), w.end(), dofIdxTarget) != w.end())
                        isVoidCol[dofIdxTarget] = true;
                }

                this->setPatternDeleteSetOfColumnsFromBCRSMatrix_(reducedMatrix, furtherReducedCoefficientMatrix[domainI][domainJ], isVoidCol);
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

     /*!
     * \brief Resizes the right-hand side
     */
    void setRHSSize(SolutionVector& RHS) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(RHS)), [&](const auto domainId)
        { RHS[domainId].resize(this->numDofs(domainId)); });
    }

    // Get the set of cell center dof indices which should be deleted from the coefficient matrix
    template<std::size_t i>
    //TODO: make it std::vector<IndexType>
    std::vector<unsigned int> reductionIndexSet(Dune::index_constant<i> domainI) const
    {
        if (i == 1)
        {
            return problem(domainI).fixedPressureScvsIndexSet();
        }
        else
        {
            // Get the set of face dof indices which should be deleted from the coefficient matri
            return problem(domainI).dirichletBoundaryScvfsIndexSet();
        }
    }

    /*!
     * \brief Resizes the reduced right-hand side
     */
    void setReducedRHSSize(SolutionVector& reducedRHS) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(reducedRHS)), [&](const auto domainId)
        { reducedRHS[domainId].resize(this->numDofs(domainId) - this->reductionIndexSet(domainId).size()); });
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
    SolutionVector& residual()
    { return *residual_; }

    //! the full reduced coefficient matrix
    JacobianMatrix& reducedCoefficientMatrix()
    { return *furtherReducedCoefficientMatrix_; }

    SolutionVector& reducedResidual()
    { return *reducedResidual_; }

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

    /*!
    * \tparam indices A set of indices that the to-be-deleted elements have in the vector.
    */
        // an example to illustrate the tought behind the loops in the following:
        // let us assume the scvfs are 0,1,...,10 with 5 and 8 being boundary scvfs
        // then we would like to delete the fifth and eightth entry of the RHS vector
        // we do so by copying new[0] = old [0], new[1] = old [1], new[2] = old [2], new[3] = old [3],
        // new[4] = old [4], new[5] = old [6], new[6] = old [7], new[7] = old [9], new[8] = old[10]
        // this could be achieved by new[i] = old [i], i=0,1,2,3,4
        // new [i-1] = old [i], i=6,7
        // new [i-2] = old [i], i=9,10
        // this means i should go over everything apart the boundary scvfs
        // and we need new [i-number of faces that we already left out]
    template<class VectorType, class IndexType>
    void removeSetOfEntriesFromVector (VectorType& vector, const std::vector<IndexType>& indices){
        if (indices.size() == 0) { return ;}

        std::vector<IndexType> tmpIndices = indices;
        std::sort (tmpIndices.begin(), tmpIndices.end());

        VectorType tmpVector;
        tmpVector.resize(vector.size() - tmpIndices.size());

        //fill intermediate reduced indices for A - delete rows
        int numBoundaryScvfsAlreadyHandled = 0;
        //k=0
        for (unsigned int i = 0; i < tmpIndices[0]; ++i){
            tmpVector[i] = vector[i];
        }
        numBoundaryScvfsAlreadyHandled ++;
        for (unsigned int k = 1; k < tmpIndices.size(); ++k){
            //k is related to the boundary scvf up to which I want to go
            for (unsigned int i = (tmpIndices[k-1]+1); i < tmpIndices[k]; ++i){
                tmpVector[i-numBoundaryScvfsAlreadyHandled] = vector[i];
            }
            numBoundaryScvfsAlreadyHandled ++;
        }
        for (unsigned int i = tmpIndices[tmpIndices.size()-1]+1; i < vector.size(); ++i){
            tmpVector[i-numBoundaryScvfsAlreadyHandled] = vector[i];
        }

        vector = tmpVector;
    }

protected:
    //! the coupling manager coupling the sub domains
    std::shared_ptr<CouplingManager> couplingManager_;

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

    // reset the reduced right-hand side vector to 0.0
    void resetReducedRHS_()
    {
        if(!reducedResidual_)
        {
            reducedResidual_ = std::make_shared<SolutionVector>();
            setReducedRHSSize(*reducedResidual_);
        }

        (*reducedResidual_) = 0.0;
    }

    // reset the reduced right-hand side vector to 0.0
    void resetReducedResidual_()
    {
        if(!reducedResidual_)
        {
            reducedResidual_ = std::make_shared<SolutionVector>();
            setReducedRHSSize(*reducedResidual_);
        }

        (*reducedResidual_) = 0.0;
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

    // reset the reduced coefficient matrix vector to 0.0
    void resetReducedCoefficientMatrix_()
    {
        if(!reducedCoefficientMatrix_ && !furtherReducedCoefficientMatrix_)
        {
            reducedCoefficientMatrix_ = std::make_shared<JacobianMatrix>();
            furtherReducedCoefficientMatrix_ = std::make_shared<JacobianMatrix>();

            setJacobianBuildMode(*reducedCoefficientMatrix_);
            setJacobianBuildMode(*furtherReducedCoefficientMatrix_);

            setReducedCoefficientMatrixPattern(*jacobian_, *reducedCoefficientMatrix_, *furtherReducedCoefficientMatrix_);
        }
        else if (!reducedCoefficientMatrix_ && furtherReducedCoefficientMatrix_)
        {
            DUNE_THROW(Dune::InvalidStateException, "Furhter reduced coefficient matrix set but reduced one not");
        }
        else if (reducedCoefficientMatrix_ && !furtherReducedCoefficientMatrix_)
        {
            DUNE_THROW(Dune::InvalidStateException, "Reduced coefficient matrix set but further reduced one not.");
        }

       (*reducedCoefficientMatrix_)  = 0.0;
       (*furtherReducedCoefficientMatrix_) = 0.0;
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
            SubDomainAssembler<i> subDomainAssembler(*this, element, curSol, *couplingManager_);
            subDomainAssembler.assembleJacobianAndResidual(jacRow, subRes, gridVariablesTuple_);
        });
    }

    template<std::size_t i, class SubRes>
    void assembleResidual_(Dune::index_constant<i> domainId, SubRes& subRes,
                           const SolutionVector& curSol)
    {
        assemble_(domainId, [&](const auto& element)
        {
            SubDomainAssembler<i> subDomainAssembler(*this, element, curSol, *couplingManager_);
            subDomainAssembler.assembleResidual(subRes);
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

    template<class MatrixType>
    void setPatternDeleteSetOfRowsFromBCRSMatrix_(const MatrixType& matrixBefore, MatrixType& matrixAfter, const std::vector<bool>& dofIsVoid){
        const std::size_t numNonVoidDofs = std::count_if(dofIsVoid.begin(), dofIsVoid.end(), [] (bool v) { return !v; });

        // reduce matrices to only dofs that take part and create index map
        std::vector<std::size_t> reductionMap(matrixBefore.N());

        std::size_t idxInReducedSpace = 0;
        for (std::size_t dofIdx = 0; dofIdx < dofIsVoid.size(); ++dofIdx)
            if (!dofIsVoid[dofIdx])
            {
                reductionMap[dofIdx] = idxInReducedSpace;
                idxInReducedSpace++;
            }

        Dune::MatrixIndexSet occupationPatternA;
        occupationPatternA.resize(numNonVoidDofs, matrixBefore.M());

        for (auto rowIt = matrixBefore.begin(); rowIt != matrixBefore.end(); ++rowIt)
            if (!dofIsVoid[rowIt.index()])
            {
                const auto reducedRowIdx = reductionMap[rowIt.index()];
                for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                    occupationPatternA.add(reducedRowIdx, colIt.index());
            }

        occupationPatternA.exportIdx(matrixAfter);
    }

    template<class MatrixType>
    void setPatternDeleteSetOfColumnsFromBCRSMatrix_(const MatrixType& matrixBefore, MatrixType& matrixAfter, const std::vector<bool>& dofIsVoid){
        const std::size_t numNonVoidDofs = std::count_if(dofIsVoid.begin(), dofIsVoid.end(), [] (bool v) { return !v; });

        // reduce matrices to only dofs that take part and create index map
        std::vector<std::size_t> reductionMap(matrixBefore.M());

        std::size_t idxInReducedSpace = 0;
        for (std::size_t dofIdx = 0; dofIdx < dofIsVoid.size(); ++dofIdx)
            if (!dofIsVoid[dofIdx])
            {
                reductionMap[dofIdx] = idxInReducedSpace;
                idxInReducedSpace++;
            }

        Dune::MatrixIndexSet occupationPatternA;
        occupationPatternA.resize(matrixBefore.N(), numNonVoidDofs);
        for (auto rowIt = matrixBefore.begin(); rowIt != matrixBefore.end(); ++rowIt)
        {
            for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
            {
                if (!dofIsVoid[colIt.index()])
                    occupationPatternA.add(rowIt.index(), reductionMap[colIt.index()]);
            }
        }

        occupationPatternA.exportIdx(matrixAfter);
    }

    //! fill the reduced coefficient matrix to 0.0
    void fillReducedCoefficientMatrix_(const JacobianMatrix& coefficientMatrix, JacobianMatrix& reducedCoefficientMatrix, JacobianMatrix& furtherReducedCoefficientMatrix)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(coefficientMatrix)), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(coefficientMatrix[domainI])), [&](const auto domainJ)
            {
                // determine the dofs that do not take part in intersections
                auto& matrix = coefficientMatrix[domainI][domainJ];
                auto& reducedMatrix = reducedCoefficientMatrix[domainI][domainJ];

                std::vector<bool> isVoidRow(matrix.N(), false);
                const auto& v=(this->reductionIndexSet(domainI));

                for (std::size_t dofIdxTarget = 0; dofIdxTarget < matrix.N(); ++dofIdxTarget)
                {
                    if (std::find(v.begin(), v.end(), dofIdxTarget) != v.end())
                        isVoidRow[dofIdxTarget] = true;
                }

                //TODO deal with empty isVoidRow or isVoidCol
                this->fillDeleteSetOfRowsFromBCRSMatrix_(matrix, reducedMatrix, isVoidRow);

                const auto& w=(this->reductionIndexSet(domainJ));
                std::vector<bool> isVoidCol(reducedMatrix.M(), false);
                for (std::size_t dofIdxTarget = 0; dofIdxTarget < reducedMatrix.M(); ++dofIdxTarget)
                {
                    if (std::find(w.begin(), w.end(), dofIdxTarget) != w.end())
                        isVoidCol[dofIdxTarget] = true;
                }

                this->fillDeleteSetOfColumnsFromBCRSMatrix_(reducedMatrix, furtherReducedCoefficientMatrix[domainI][domainJ], isVoidCol);
            });
        });
    }

    template<class MatrixType>
    void fillDeleteSetOfRowsFromBCRSMatrix_(const MatrixType& matrixBefore, MatrixType& matrixAfter, const std::vector<bool>& dofIsVoid)
    {
        const std::size_t numNonVoidDofs = std::count_if(dofIsVoid.begin(), dofIsVoid.end(), [] (bool v) { return !v; });

        // reduce matrices to only dofs that take part and create index map
        std::vector<std::size_t> reductionMap(matrixBefore.N());

        std::size_t idxInReducedSpace = 0;
        for (std::size_t dofIdx = 0; dofIdx < dofIsVoid.size(); ++dofIdx)
            if (!dofIsVoid[dofIdx])
            {
                reductionMap[dofIdx] = idxInReducedSpace;
                idxInReducedSpace++;
            }

        for (auto rowIt = matrixBefore.begin(); rowIt != matrixBefore.end(); ++rowIt)
            if (!dofIsVoid[rowIt.index()])
            {
                const auto reducedRowIdx = reductionMap[rowIt.index()];
                for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                    matrixAfter[reducedRowIdx][colIt.index()] = *colIt;
            }
    }

    template<class MatrixType>
    void fillDeleteSetOfColumnsFromBCRSMatrix_(const MatrixType& matrixBefore, MatrixType& matrixAfter, const std::vector<bool>& dofIsVoid){
        const std::size_t numNonVoidDofs = std::count_if(dofIsVoid.begin(), dofIsVoid.end(), [] (bool v) { return !v; });

        // reduce matrices to only dofs that take part and create index map
        std::vector<std::size_t> reductionMap(matrixBefore.M());

        std::size_t idxInReducedSpace = 0;
        for (std::size_t dofIdx = 0; dofIdx < dofIsVoid.size(); ++dofIdx)
            if (!dofIsVoid[dofIdx])
            {
                reductionMap[dofIdx] = idxInReducedSpace;
                idxInReducedSpace++;
            }

        for (auto rowIt = matrixBefore.begin(); rowIt != matrixBefore.end(); ++rowIt)
        {
            for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                if (!dofIsVoid[colIt.index()])
                    matrixAfter[rowIt.index()][reductionMap[colIt.index()]] = *colIt;
        }
    }

    void fillReducedRHS_(const SolutionVector& RHS, SolutionVector& reducedRHS)
    {
        reducedRHS = RHS;
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(reducedRHS)), [&](const auto domainId)
        { this->removeSetOfEntriesFromVector(reducedRHS[domainId], this->reductionIndexSet(domainId)); });
    }

    void fillReducedResidual_(const SolutionVector& residual, SolutionVector& reducedResidual)
    {
        reducedResidual = residual;
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(reducedResidual)), [&](const auto domainId)
        { this->removeSetOfEntriesFromVector(reducedResidual[domainId], this->reductionIndexSet(domainId)); });
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
    std::shared_ptr<SolutionVector> residual_;
    std::shared_ptr<JacobianMatrix> reducedCoefficientMatrix_;
    std::shared_ptr<SolutionVector> reducedResidual_;
    std::shared_ptr<JacobianMatrix> furtherReducedCoefficientMatrix_;
};

} // end namespace Dumux

#endif
