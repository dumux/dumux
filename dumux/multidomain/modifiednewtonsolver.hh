// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Nonlinear
 * \brief Modified Newton Solver.
 */
#ifndef DUMUX_MULTIDOMAIN_MODIFIED_NEWTON_SOLVER_HH
#define DUMUX_MULTIDOMAIN_MODIFIED_NEWTON_SOLVER_HH

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/porenetwork/2p/newtonconsistencychecks.hh>

namespace Dumux {
namespace Detail {

template <typename T>
using InvasionStateDetector = decltype(std::declval<T>().invasionState());

template<class T>
static constexpr bool hasInvasionState()
{ return Dune::Std::is_detected<InvasionStateDetector, T>::value; }


// Primary template
template <class VolumeVariables, typename U = int>
struct IsMPNC : std::false_type
{};

// Specialization for MPNC
template <class VolumeVariables>
struct IsMPNC <VolumeVariables, decltype((void) VolumeVariables::Indices::s0Idx, 0)> : std::true_type
{};

// Primary template
template <class VolumeVariables, typename U = int>
struct IsNonIsothermal : std::false_type
{};

// Specialization for non-isothermal models
template <class VolumeVariables>
struct IsNonIsothermal <VolumeVariables, decltype((void) VolumeVariables::Indices::temperatureIdx, 0)> : std::true_type
{};

//! This is a partial copy of Dune::MultiTypeBlockMatrix. TODO: It can be deleted once Dune::MultiTypeBlockMatrix
//! exposes std::tuple's constructor.
template<class FirstRow, class... Args>
class MultiTypeBlockMatrix : std::tuple<FirstRow, Args...>
{
    using ParentType = std::tuple<FirstRow, Args...>;
  public:

    using ParentType::ParentType;

    /**
     * own class' type
     */
    using type = MultiTypeBlockMatrix<FirstRow, Args...>;

    /** \brief Type used for sizes */
    using size_type = std::size_t;

    using field_type = typename FirstRow::field_type;

    /** \brief Return the number of matrix rows */
    static constexpr size_type N()
    {
      return 1+sizeof...(Args);
    }

    /** \brief Return the number of matrix columns */
    static constexpr size_type M()
    {
      return FirstRow::size();
    }

    template< size_type index >
    auto
    operator[] ( const std::integral_constant< size_type, index > indexVariable ) -> decltype(std::get<index>(*this))
    {
      DUNE_UNUSED_PARAMETER(indexVariable);
      return std::get<index>(*this);
    }

   /** \brief Const random-access operator
     *
     * This is the const version of the random-access operator.  See the non-const version for a full
     * explanation of how to use it.
     */
    template< size_type index >
    auto
    operator[] ( const std::integral_constant< size_type, index > indexVariable ) const -> decltype(std::get<index>(*this))
    {
      DUNE_UNUSED_PARAMETER(indexVariable);
      return std::get<index>(*this);
    }

    /** \brief y = A x
     */
    template<typename X, typename Y>
    void mv (const X& x, Y& y) const {
      static_assert(X::size() == M(), "length of x does not match row length");
      static_assert(Y::size() == N(), "length of y does not match row count");
      y = 0; //reset y (for mv uses umv)
      umv(x,y);
    }

    /** \brief y += A x
     */
    template<typename X, typename Y>
    void umv (const X& x, Y& y) const {
      static_assert(X::size() == M(), "length of x does not match row length");
      static_assert(Y::size() == N(), "length of y does not match row count");
      using namespace Dune::Hybrid;
      forEach(integralRange(Dune::Hybrid::size(y)), [&](auto&& i) {
        using namespace Dune::Hybrid; // needed for icc, see issue #31
        forEach(integralRange(Dune::Hybrid::size(x)), [&](auto&& j) {
          (*this)[i][j].umv(x[j], y[i]);
        });
      });
    }

    /** \brief y += alpha A x
     */
    template<typename AlphaType, typename X, typename Y>
    void usmv (const AlphaType& alpha, const X& x, Y& y) const {
      static_assert(X::size() == M(), "length of x does not match row length");
      static_assert(Y::size() == N(), "length of y does not match row count");
      using namespace Dune::Hybrid;
      forEach(integralRange(Dune::Hybrid::size(y)), [&](auto&& i) {
        using namespace Dune::Hybrid; // needed for icc, see issue #31
        forEach(integralRange(Dune::Hybrid::size(x)), [&](auto&& j) {
          (*this)[i][j].usmv(alpha, x[j], y[i]);
        });
      });
    }

};

/*!
 * \brief a function to get a MultiTypeBlockVector with const references to some entries of another MultiTypeBlockVector
 * \param v a MultiTypeBlockVector
 * \param indices the indices of the entries that should be referenced
 * TODO can be removed if it gets implemented in dumux-master
 */
template<class ...Args, std::size_t ...i>
auto partial(const Dune::MultiTypeBlockVector<Args...>& v, Dune::index_constant<i>... indices)
{
    return Dune::MultiTypeBlockVector<std::add_lvalue_reference_t<const std::decay_t<std::tuple_element_t<indices, std::tuple<Args...>>>>...>(v[indices]...);
}

} // end namespace Detail

/*!
 * \ingroup Nonlinear
 * \brief An implementation of a Newton solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note If you want to specialize only some methods but are happy with the
 *       defaults of the reference solver, derive your solver from
 *       this class and simply overload the required methods.
 */
template <class Assembler, class LinearSolver, class CouplingManager,
          class Reassembler = DefaultPartialReassembler,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>,
          template<class, class> class NewtonConsistencyChecks = PoreNetwork::TwoPNewtonConsistencyChecks>
class ModifiedMultiDomainNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>
{
    using ParentType = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>;
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SolutionVector = typename Assembler::ResidualType;

public:
    using ParentType::ParentType;

    /*!
     * \brief Called after each Newton update
     *
     * \param uCurrentIter The current global solution vector
     * \param uLastIter The previous global solution vector
     */
     void newtonEndStep(SolutionVector &uCurrentIter,
                        const SolutionVector &uLastIter) final
    {
        // call the method of the base class
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        ++numSteps_;

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& gridVariables = this->assembler().gridVariables(id);
            const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

            updateInvsionState_(uCurrentIter, id, std::integral_constant<bool, hasInvasionState>());
        });
    }

    /*!
     * \brief Returns true if the current solution can be considered to
     *        be accurate enough
     */
    bool newtonConverged() const final
    {
        if (Dune::any_true(switchedInLastIteration_))
            return false;

        if (lastIterationWasChopped_)
            return false;

        return ParentType::newtonConverged();
    }

    void newtonFail(SolutionVector& u) final
    {
        ParentType::newtonFail(u);

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& gridVariables = this->assembler().gridVariables(id);
            const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

            resetInvasionState_(id, std::integral_constant<bool, hasInvasionState>());
        });

        if (invasionEventInTimeStep_)
        {
            std::cout << "Newton failed due to invasion event" << std::endl;
            static const Scalar invasionFactor = getParam<Scalar>("Newton.InvasionEventReductionFactor", 1.0);
            this->setRetryTimeStepReductionFactor(origRetryReductionFactor_ * invasionFactor);
        }
        else
            this->setRetryTimeStepReductionFactor(origRetryReductionFactor_);

        invasionEventInTimeStep_ = false;
        lastIterationWasChopped_ = false;
    }

    /*!
     * \brief Called if the Newton method ended successfully
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed() final
    {
        ParentType::newtonSucceed();

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& gridVariables = this->assembler().gridVariables(id);
            const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();

            advanceInvasionState_(id, std::integral_constant<bool, hasInvasionState>());
        });
    }

    /*!
     *
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    void newtonBegin(SolutionVector& u) final
    {
        ParentType::newtonBegin(u);
        invasionEventInTimeStep_ = false;
        numSteps_ = 0;
        lastIterationWasChopped_ = false;
    }


private:

    /*!
     * \brief Update the current solution of the Newton method
     *
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void choppedUpdate_(SolutionVector &uCurrentIter,
                        const SolutionVector &uLastIter,
                        const SolutionVector &deltaU) final
    {
        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;

        static const int maxNumChoppedUpdates = getParam<int>("Newton.NumChoppedUpdates", 0);

        if (maxNumChoppedUpdates > 0 && numSteps_ <= maxNumChoppedUpdates * (invasionEventInTimeStep_ ? 2 : 1))
        {
            lastIterationWasChopped_ = true;

            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](const auto id)
            {
                auto& gridVariables = this->assembler().gridVariables(id);
                const auto hasInvasionState = Detail::hasInvasionState<std::decay_t< decltype(gridVariables.gridFluxVarsCache())>>();
                using VolumeVariables = typename std::decay_t<decltype(gridVariables)>::GridVolumeVariables::VolumeVariables;
                using Indices = typename VolumeVariables::Indices;



                if constexpr (hasInvasionState)
                {
                    static const bool useTransFactor = getParam<bool>("Newton.UseTransFactor", false);
                    if (useTransFactor && invasionEventInTimeStep_)
                    {
                        static const Scalar firstFactor = getParam<Scalar>("Newton.FirstTransFactor", 1e-5);
                        static const Scalar progression = getParam<Scalar>("Newton.TransfactorProgression", 1.0);
                        gridVariables.gridFluxVarsCache().invasionState().setTransmissibilityFactor(firstFactor*(numSteps_*progression + 1));
                        std::cout << "Factor is " << gridVariables.gridFluxVarsCache().invasionState().transmissibilityFactor();
                    }

                    if constexpr (Detail::IsMPNC<VolumeVariables>::value)
                    {
                        std::cout << "using MPNC chopped update" << std::endl;
                        static constexpr auto numPhases = VolumeVariables::numFluidPhases();
                        static constexpr auto numComponents = VolumeVariables::numFluidComponents();

                        for (std::size_t i = 0; i < uCurrentIter[id].size(); ++i)
                        {
                            auto& uCurr = uCurrentIter[id][i];
                            const auto& uLast = uLastIter[id][i];
                            for (std::size_t phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                                saturationChop_(uCurr[Indices::s0Idx + phaseIdx],
                                                uLast[Indices::s0Idx + phaseIdx]);

                            pressureChop_(uCurr[Indices::p0Idx],
                                          uLast[Indices::p0Idx]);

                            for (std::size_t comp = 0; comp < numComponents; ++comp)
                                pressureChop_(uCurr[Indices::fug0Idx + comp],
                                              uLast[Indices::fug0Idx + comp]);

                            if constexpr (Detail::IsNonIsothermal<VolumeVariables>::value)
                                temperatureChop_(uCurr[Indices::temperatureIdx],
                                                 uLast[Indices::temperatureIdx]);
                        }
                    }
                    else
                    {
                        std::cout << "using 2pnc chopped update" << std::endl;
                        for (std::size_t i = 0; i < uCurrentIter[id].size(); ++i)
                        {
                            auto& uCurr = uCurrentIter[id][i];
                            const auto& uLast = uLastIter[id][i];

                            saturationChop_(uCurr[Indices::switchIdx],
                                            uLast[Indices::switchIdx]);
                            pressureChop_(uCurr[Indices::pressureIdx],
                                          uLast[Indices::pressureIdx]);

                            if constexpr (Detail::IsNonIsothermal<VolumeVariables>::value)
                                temperatureChop_(uCurr[Indices::temperatureIdx],
                                                 uLast[Indices::temperatureIdx]);
                        }
                    }
                }
//                 else if constexpr (id == 1)
//                 {
//                     for (std::size_t i = 0; i < uCurrentIter[id].size(); ++i)
//                     {
//                         // chop the free flow, too
//                         pressureChop_(uCurrentIter[id][i][Indices::pressureIdx],
//                                       uLastIter[id][i][Indices::pressureIdx]);
//                     }
//                 }
            });
        }
        else
            lastIterationWasChopped_ = false;


        if (this->enableResidualCriterion())
            this->computeResidualReduction_(uCurrentIter);

        else
        {
            // If we get here, the convergence criterion does not require
            // additional residual evalutions. Thus, the grid variables have
            // not yet been updated to the new uCurrentIter.
            this->assembler().updateGridVariables(uCurrentIter);
        }
    }

    static void clampValue_(Scalar &val,
                            const Scalar minVal,
                            const Scalar maxVal)
    {
        using std::max;
        using std::min;
        val = max(minVal, min(val, maxVal));
    };

    static void pressureChop_(Scalar &val,
                              const Scalar oldVal)
    {
        using std::max; using std::clamp;
        static const Scalar maxDeltaInput = getParam<Scalar>("Newton.PressureMaxDelta", -1.0);
        const Scalar maxDelta = maxDeltaInput > 0.0 ? maxDeltaInput : max(oldVal/4.0, 10e3);
        val = clamp(val, oldVal - maxDelta, oldVal + maxDelta);
        val = max(0.0, val); // don't allow negative pressures
    }

    static void saturationChop_(Scalar &val,
                                const Scalar oldVal)
    {
        static const Scalar maxDelta = getParam<Scalar>("Newton.SaturationMaxDelta", 0.25);
        using std::clamp;
        val = clamp(val, oldVal - maxDelta, oldVal + maxDelta);
        static const auto range = getParam<std::array<Scalar, 2>>("Newton.SaturationClampValues", std::array<Scalar, 2>{0.0, 1.0});
        val = clamp(val, range[0], range[1]);
    }

    static void temperatureChop_(Scalar &val,
                                 const Scalar oldVal)
    {
        static const Scalar maxDelta = getParam<Scalar>("Newton.TemperatureMaxDelta", 2.0);
        using std::clamp;
        val = clamp(val, oldVal - maxDelta, oldVal + maxDelta);
        static const auto range = getParam<std::array<Scalar, 2>>("Newton.TemperatureClampValues", std::array<Scalar, 2>{273.0, 600.00});
        val = clamp(val, range[0], range[1]);
    }


    template<std::size_t i>
    void updateInvsionState_(const SolutionVector& uCurrentIter, Dune::index_constant<i> id, std::true_type)
    {
        std::cout << "updating invasion state for " << i << std::endl;
        auto& gridVariables = this->assembler().gridVariables(id);
        switchedInLastIteration_[i] =  gridVariables.gridFluxVarsCache().invasionState().update(uCurrentIter[id],
                                                                                                gridVariables.curGridVolVars(),
                                                                                                gridVariables.gridFluxVarsCache());

        if (switchedInLastIteration_[i])
        {
            invasionEventInTimeStep_ = true;
            numSteps_ = 0;
        }

        // If the solution is about to be accepted, check for accuracy and trigger a retry
        // with a decreased time step size if necessary.
        if (newtonConverged())
        {
            PoreNetwork::TwoPNewtonConsistencyChecks<std::decay_t<decltype(gridVariables)>, std::decay_t<decltype(uCurrentIter[id])>> checks;
            checks.performChecks(gridVariables, uCurrentIter[id], this->assembler().prevSol()[id]);
        }
    }

    template<std::size_t i>
    void resetInvasionState_(Dune::index_constant<i> id, std::true_type)
    {
        this->assembler().gridVariables(id).gridFluxVarsCache().invasionState().reset();
    }

    template<std::size_t i>
    void advanceInvasionState_(Dune::index_constant<i> id, std::true_type)
    {
        this->assembler().gridVariables(id).gridFluxVarsCache().invasionState().advance();
    }

    template<std::size_t i>
    void updateInvsionState_(const SolutionVector& uCurrentIter, Dune::index_constant<i> id, std::false_type)
    {
        switchedInLastIteration_[i] = false;
    }

    template<std::size_t i>
    void resetInvasionState_(Dune::index_constant<i> id, std::false_type)
    {}

    template<std::size_t i>
    void advanceInvasionState_(Dune::index_constant<i> id, std::false_type)
    {}


    bool solveLinearSystem_(SolutionVector& deltaU) final
    {
        auto& ls = this->linearSolver();
        auto& A = this->assembler().jacobian();
        auto& x = deltaU;
        auto& b = this->assembler().residual();

        assert(this->checkSizesOfSubMatrices(A) && "Sub-blocks of MultiTypeBlockMatrix have wrong sizes!");

        using namespace Dune::Indices;

        const auto& A00 = A[_0][_0];
        const auto& A11 = A[_1][_1];
        const auto& A22 = A[_2][_2];

        const auto& A01 = A[_0][_1];
        const auto& A02 = A[_0][_2];

        const auto& A10 = A[_1][_0];
        const auto& A12 = A[_1][_2];

        const auto& A20 = A[_2][_0];
        const auto& A21 = A[_2][_1];


        using FirstRow = Dune::MultiTypeBlockVector<const std::decay_t<decltype(A11)>&, const std::decay_t<decltype(A10)>&, const std::decay_t<decltype(A12)>&>;
        using SecondRow = Dune::MultiTypeBlockVector<const std::decay_t<decltype(A01)>&, const std::decay_t<decltype(A00)>&, const std::decay_t<decltype(A02)>&>;
        using ThirdRow = Dune::MultiTypeBlockVector<const std::decay_t<decltype(A21)>&, const std::decay_t<decltype(A20)>&, const std::decay_t<decltype(A22)>&>;

        using NewA = Detail::MultiTypeBlockMatrix<FirstRow, SecondRow, ThirdRow>;

        FirstRow firstRow(A11, A10, A12);
        SecondRow secondRow(A01, A00, A02);
        ThirdRow thirdRow(A21, A20, A22);

        const NewA newA = NewA(firstRow, secondRow, thirdRow);

        const auto newM = MatrixConverter<NewA>::multiTypeToBCRSMatrix(newA);

        using namespace Dune::Indices;

        auto newB = partial(b, _1, _0, _2);
        const auto newBtmp = VectorConverter<decltype(newB)>::multiTypeToBlockVector(newB);
        auto newX = partial(x, _1, _0, _2);
        const std::size_t numRows = newM.N();

        // create a blockvector to which the linear solver writes the solution
        using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
        using BlockVector = typename Dune::BlockVector<VectorBlock>;
        BlockVector y(numRows);

        // solve
        const bool converged = ls.solve(newM, y, newBtmp);

        // copy back the result y into x
        if(converged)
            VectorConverter<decltype(newB)>::retrieveValues(newX, y);

        return converged;
    }

    std::array<bool, Assembler::Traits::numSubDomains> switchedInLastIteration_ = {};
    Scalar origRetryReductionFactor_ = this->retryTimeStepReductionFactor();
    bool invasionEventInTimeStep_ = false;
    int numSteps_ = 0;
    bool lastIterationWasChopped_ = false;
};

} // end namespace Dumux

#endif
