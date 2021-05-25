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
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 */
#ifndef DUMUX_FV_ASSEMBLER_HH
#define DUMUX_FV_ASSEMBLER_HH

#include "tbb/parallel_for.h"

#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/io/format.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/method.hh>
#include <dumux/linear/parallelhelpers.hh>

#include "jacobianpattern.hh"
#include "diffmethod.hh"
#include "boxlocalassembler.hh"
#include "cclocalassembler.hh"
#include "fclocalassembler.hh"

namespace Dumux::Detail {

template<DiscretizationMethod diffMethod>
struct LocalAssemblerChooser;

template<>
struct LocalAssemblerChooser<DiscretizationMethod::box>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod, bool isImplicit>
    using type = BoxLocalAssembler<TypeTag, Impl, diffMethod, isImplicit>;
};

template<>
struct LocalAssemblerChooser<DiscretizationMethod::ccmpfa>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod, bool isImplicit>
    using type = CCLocalAssembler<TypeTag, Impl, diffMethod, isImplicit>;
};

template<>
struct LocalAssemblerChooser<DiscretizationMethod::cctpfa>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod, bool isImplicit>
    using type = CCLocalAssembler<TypeTag, Impl, diffMethod, isImplicit>;
};

template<>
struct LocalAssemblerChooser<DiscretizationMethod::fcstaggered>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod, bool isImplicit>
    using type = FaceCenteredLocalAssembler<TypeTag, Impl, diffMethod, isImplicit>;
};

template<class TypeTag, class Impl, DiffMethod diffMethod, bool isImplicit>
using LocalAssemblerChooser_t = typename LocalAssemblerChooser<
    GetPropType<TypeTag, Properties::GridGeometry>::discMethod
>::template type<TypeTag, Impl, diffMethod, isImplicit>;

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam isImplicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
class FVAssembler
{
    using GridGeo = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeo::GridView;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename GridView::Grid::template Codim<0>::EntitySeed;
    using TimeLoop = TimeLoopBase<GetPropType<TypeTag, Properties::Scalar>>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr bool isBox = GridGeo::discMethod == DiscretizationMethod::box;

    using ThisType = FVAssembler<TypeTag, diffMethod, isImplicit>;
    using LocalAssembler = typename Detail::LocalAssemblerChooser_t<TypeTag, ThisType, diffMethod, isImplicit>;

public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridGeometry = GridGeo;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using ResidualType = SolutionVector;

    /*!
     * \brief The constructor for stationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    FVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , timeLoop_()
    , isStationaryProblem_(true)
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
    }

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    FVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<GridVariables> gridVariables,
                std::shared_ptr<const TimeLoop> timeLoop,
                const SolutionVector& prevSol)
    : problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , timeLoop_(timeLoop)
    , prevSol_(&prevSol)
    , isStationaryProblem_(!timeLoop)
    {}

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const SolutionVector& curSol, const PartialReassembler* partialReassembler = nullptr)
    {
        checkAssemblerState_();
        resetJacobian_(partialReassembler);
        resetResidual_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, *gridVariables_, partialReassembler);
        });

        enforcePeriodicConstraints_(*jacobian_, *residual_, curSol, *gridGeometry_);
    }

    /*!
     * \brief Assembles only the global Jacobian of the residual.
     */
    void assembleJacobian(const SolutionVector& curSol)
    {
        checkAssemblerState_();
        resetJacobian_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *gridVariables_);
        });
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        assembleResidual(*residual_, curSol);
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const SolutionVector& curSol) const
    {
        checkAssemblerState_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleResidual(r);
        });
    }

    //! compute the residual and return it's vector norm
    Scalar residualNorm(const SolutionVector& curSol) const
    {
        ResidualType residual(numDofs());
        assembleResidual(residual, curSol);

        // issue a warning if the caluclation is used in parallel with overlap
        static bool warningIssued = false;

        if (gridView().comm().size() > 1 && gridView().overlapSize(0) == 0)
        {
            if constexpr (isBox)
            {
                using DM = typename GridGeometry::VertexMapper;
                using PVHelper = ParallelVectorHelper<GridView, DM, GridView::dimension>;

                PVHelper vectorHelper(gridView(), gridGeometry_->vertexMapper());

                vectorHelper.makeNonOverlappingConsistent(residual);
            }
        }
        else if (!warningIssued)
        {
            if (gridView().comm().size() > 1 && gridView().comm().rank() == 0)
                std::cout << "\nWarning: norm calculation adds entries corresponding to\n"
                << "overlapping entities multiple times. Please use the norm\n"
                << "function provided by a linear solver instead." << std::endl;

            warningIssued = true;
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

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

        // check and/or set the BCRS matrix's build mode
        if (jacobian_->buildMode() == JacobianMatrix::BuildMode::unknown)
            jacobian_->setBuildMode(JacobianMatrix::random);
        else if (jacobian_->buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        setJacobianPattern();
        setResidualSize();
        buildElementSets_();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        jacobian_->setBuildMode(JacobianMatrix::random);
        residual_ = std::make_shared<SolutionVector>();

        setJacobianPattern();
        setResidualSize();
        buildElementSets_();
    }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // create occupation pattern of the jacobian
        const auto occupationPattern = getJacobianPattern<isImplicit>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return gridGeometry_->numDofs(); }

    //! The problem
    const Problem& problem() const
    { return *problem_; }

    //! The global finite volume geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The gridview
    const GridView& gridView() const
    { return gridGeometry().gridView(); }

    //! The global grid variables
    GridVariables& gridVariables()
    { return *gridVariables_; }

    //! The global grid variables
    const GridVariables& gridVariables() const
    { return *gridVariables_; }

    //! The jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! The residual vector (rhs)
    SolutionVector& residual()
    { return *residual_; }

    //! The solution of the previous time step
    const SolutionVector& prevSol() const
    { return *prevSol_; }

    /*!
     * \brief Set time loop for instationary problems
     * \note calling this turns this into a stationary assembler
     */
    void setTimeLoop(std::shared_ptr<const TimeLoop> timeLoop)
    { timeLoop_ = timeLoop; isStationaryProblem_ = !static_cast<bool>(timeLoop); }

    /*!
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { prevSol_ = &u;  }

    /*!
     * \brief Whether we are assembling a stationary or instationary problem
     */
    bool isStationaryProblem() const
    { return isStationaryProblem_; }

    /*!
     * \brief Create a local residual object (used by the local assembler)
     */
    LocalResidual localResidual() const
    { return LocalResidual(problem_.get(), timeLoop_.get()); }

    /*!
     * \brief Update the grid variables
     */
    void updateGridVariables(const SolutionVector &cursol)
    {
        gridVariables().update(cursol);
    }

    /*!
     * \brief Reset the gridVariables
     */
    void resetTimeStep(const SolutionVector &cursol)
    {
        gridVariables().resetTimeStep(cursol);
    }

private:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize();
        }

        (*residual_) = 0.0;
    }

    // reset the Jacobian matrix to 0.0
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        if (partialReassembler)
            partialReassembler->resetJacobian(*this);
        else
            *jacobian_ = 0.0;
    }

    // check if the assembler is in a correct state for assembly
    void checkAssemblerState_() const
    {
        if (!isStationaryProblem_ && !prevSol_)
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     */
    template<typename AssembleElementFunc>
    void assemble_(AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
        try
        {
            // make this element loop run in parallel
            // for this we have to color the elements so that we don't get
            // race conditions when writing into the global matrix
            // each color can be assembled using multiple threads
            // this is because locking with mutexes is expensive and
            // we assume that we still have enough elements per color
            // to make use of multiple cores

            // let the local assembler add the element contributions
            if (parallelAssembly_)
            {
                for (const auto& elements : elementSets_)
                {
                    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, elements.size()),
                    [&](const tbb::blocked_range<size_t>& r)
                    {
                        for (std::size_t i=r.begin(); i<r.end(); ++i)
                        {
                            const auto element = gridView().grid().entity(elements[i]);
                            assembleElement(element);
                        }
                    });
                }
            }
            else
            {
                for (const auto& element : elements(gridView()))
                    assembleElement(element);
            }

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView().comm().size() > 1)
            succeeded = gridView().comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    template<class GG> std::enable_if_t<GG::discMethod == DiscretizationMethod::box, void>
    enforcePeriodicConstraints_(JacobianMatrix& jac, SolutionVector& res, const SolutionVector& curSol, const GG& gridGeometry)
    {
        for (const auto& m : gridGeometry.periodicVertexMap())
        {
            if (m.first < m.second)
            {
                // add the second row to the first
                res[m.first] += res[m.second];
                const auto end = jac[m.second].end();
                for (auto it = jac[m.second].begin(); it != end; ++it)
                    jac[m.first][it.index()] += (*it);

                // enforce constraint in second row
                res[m.second] = curSol[m.second] - curSol[m.first];
                for (auto it = jac[m.second].begin(); it != end; ++it)
                    (*it) = it.index() == m.second ? 1.0 : it.index() == m.first ? -1.0 : 0.0;
            }
        }
    }

    template<class GG> std::enable_if_t<GG::discMethod != DiscretizationMethod::box, void>
    enforcePeriodicConstraints_(JacobianMatrix& jac, SolutionVector& res, const SolutionVector& curSol, const GG& gridGeometry) {}

    void buildElementSets_()
    {
        parallelAssembly_ = getParam<bool>("Assembly.Parallel", false);
        if (parallelAssembly_)
        {
            Dune::Timer timer;
            // cell-centered algorithm
            // greedy algorithm
            std::vector<int> localColors; localColors.reserve(10);
            std::vector<bool> notAssigned; notAssigned.reserve(10);
            std::vector<int> colors(gridView().size(0), -1);
            for (const auto& element : elements(gridView()))
            {
                localColors.clear();
                for (const auto& i : intersections(gridView(), element))
                {
                    if (i.neighbor())
                    {
                        const auto& neighbor = i.outside();
                        localColors.push_back(colors[gridGeometry().elementMapper().index(neighbor)]);
                        for (const auto& j : intersections(gridView(), neighbor))
                            if (j.neighbor())
                                localColors.push_back(colors[gridGeometry().elementMapper().index(j.outside())]);
                    }
                }

                // find smallest color (positive integer) not in localColors
                const auto smallestAvailableColor = [&]
                {
                    const int numLocalColors = localColors.size();
                    notAssigned.assign(numLocalColors, true);

                    // worst case for numLocalColors=3 is localColors={0, 1, 2}
                    // which should result in 3 as smallest available color
                    for (int i = 0; i < numLocalColors; i++)
                        if (localColors[i] >= 0 && localColors[i] < numLocalColors)
                            notAssigned[localColors[i]] = false;

                    for (int i = 0; i < numLocalColors; i++)
                        if (notAssigned[i])
                            return i;

                    return numLocalColors;
                }();

                colors[gridGeometry().elementMapper().index(element)]
                    = smallestAvailableColor;

                if (smallestAvailableColor < elementSets_.size())
                    elementSets_[smallestAvailableColor].push_back(element.seed());
                else
                    elementSets_.push_back(std::vector<ElementSeed>{ element.seed() });
            }

            std::cout << Fmt::format("Colored elements with {} colors in {} seconds.\n",
                                     elementSets_.size(), timer.elapsed());
        }
    }

    //! pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    //! the finite volume geometry of the grid
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    //! the time loop for instationary problem assembly
    std::shared_ptr<const TimeLoop> timeLoop_;

    //! an observing pointer to the previous solution for instationary problems
    const SolutionVector* prevSol_ = nullptr;

    //! if this assembler is assembling an instationary problem
    bool isStationaryProblem_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<SolutionVector> residual_;

    //! element sets for parallel assembly
    bool parallelAssembly_ = false;
    std::deque<std::vector<ElementSeed>> elementSets_;
};

} // namespace Dumux

#endif
