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
 * \brief An assembler for Jacobians and residuals of PDEs defined on grids.
 */
#ifndef DUMUX_ASSEMBLY_ASSEMBLER_HH
#define DUMUX_ASSEMBLY_ASSEMBLER_HH

#include <memory>
#include <utility>
#include <concepts>

#include <dumux/experimental/new_assembly/dumux/assembly/operatorweights.hh>
#include <dumux/experimental/new_assembly/dumux/common/variables.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename T>
concept LocalAssembler = requires {
    typename T::GridGeometry;
    typename T::GridVariables;
    Concepts::Variables<typename T::GridVariables>;
};

template<typename T>
concept BasicLocalAssembler
    = LocalAssembler<T>
    and std::constructible_from<T,
                                const typename T::GridGeometry::LocalView&,
                                typename T::GridVariables::LocalView&>;

template<typename T>
concept WeightedOperatorsLocalAssembler
    = LocalAssembler<T>
    and std::constructible_from<T,
                                const typename T::GridGeometry::LocalView&,
                                typename T::GridVariables::LocalView&,
                                const OperatorWeights<Variables::ScalarType<typename T::GridVariables>>&>;

template<typename LA, typename J>
concept LocalJacobianAssembler = requires(LA& localAssembler, J& jacobian) {
    { localAssembler.addJacobianEntries(jacobian) };
};

template<typename LA, typename R>
concept LocalResidualAssembler = requires(LA& localAssembler, R& residual) {
    { localAssembler.addResidualEntries(residual) };
};

template<typename LA, typename J, typename R>
concept LocalJacobianAndResidualAssembler = requires(LA& localAssembler,
                                                     J& jacobian,
                                                     R& residual) {
    { localAssembler.addJacobianAndResidualEntries(jacobian, residual) };
};

} // namespace Detail
#endif // DOXYGEN


/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobians and residuals of PDEs defined on grids.
 */
template<typename LA>
class Assembler
{
    using GridVariables = typename LA::GridVariables;
    using Scalar = Variables::ScalarType<GridVariables>;

public:
    using LocalAssembler = LA;
    using Variables = GridVariables;
    using GridGeometry = typename LA::GridGeometry;
    using OperatorWeights = Dumux::OperatorWeights<Scalar>;

    explicit Assembler(std::shared_ptr<const GridGeometry> gg)
    : gridGeometry_(std::move(gg))
    {}

    //! Jacobian assembly
    template<typename JacobianMatrix> requires(
        Detail::BasicLocalAssembler<LA> and
        Detail::LocalJacobianAssembler<LA, JacobianMatrix>)
    void assembleJacobian(JacobianMatrix& jacobian, const Variables& variables) const
    {
        assemble_(variables, [&] (const auto& localGeometry, auto& localVariables) {
            LocalAssembler{localGeometry, localVariables}.addJacobianEntries(jacobian);
        });
    }

    //! Residual assembly
    template<typename Residual> requires(
        Detail::BasicLocalAssembler<LA> and
        Detail::LocalResidualAssembler<LA, Residual>)
    void assembleResidual(Residual& residual, const Variables& variables) const
    {
        assemble_(variables, [&] (const auto& localGeometry, auto& localVariables) {
            LA{localGeometry, localVariables}.addResidualEntries(residual);
        });
    }

    //! Jacobian & residual assembly
    template<typename JacobianMatrix, typename Residual> requires(
        Detail::BasicLocalAssembler<LA> and
        Detail::LocalJacobianAndResidualAssembler<LA, JacobianMatrix, Residual>)
    void assembleJacobianAndResidual(JacobianMatrix& jacobian,
                                     Residual& residual,
                                     const Variables& variables) const
    {
        assemble_(variables, [&] (const auto& localGeometry, auto& localVariables) {
            LA{localGeometry, localVariables}.addJacobianAndResidualEntries(jacobian, residual);
        });
    }

    //! Jacobian assembly with weights for the spatial/temporal operators
    template<typename JacobianMatrix> requires(
        Detail::WeightedOperatorsLocalAssembler<LA> and
        Detail::LocalJacobianAssembler<LA, JacobianMatrix>)
    void assembleJacobian(JacobianMatrix& jacobian,
                          const Variables& variables,
                          const OperatorWeights& weights) const
    {
        assemble_(variables, [&] (const auto& localGeometry, auto& localVariables) {
            LA{localGeometry, localVariables, weights}.addJacobianEntries(jacobian);
        });
    }

    //! Residual assembly with weights for temporal/spatial operators
    template<typename Residual> requires(
        Detail::WeightedOperatorsLocalAssembler<LA> and
        Detail::LocalResidualAssembler<LA, Residual>)
    void assembleResidual(Residual& residual,
                          const Variables& variables,
                          const OperatorWeights& weights) const
    {
        assemble_(variables, [&] (const auto& localGeometry, auto& localVariables) {
            LA{localGeometry, localVariables, weights}.addResidualEntries(residual);
        });
    }

    //! Jacobian & residual assembly with weights for temporal/spatial operators
    template<typename JacobianMatrix, typename Residual> requires(
        Detail::WeightedOperatorsLocalAssembler<LA> and
        Detail::LocalJacobianAndResidualAssembler<LA, JacobianMatrix, Residual>)
    void assembleJacobianAndResidual(JacobianMatrix& jacobian,
                                     Residual& residual,
                                     const Variables& variables,
                                     const OperatorWeights& weights) const
    {
        assemble_(variables, [&] (const auto& localGeometry, auto& localVariables) {
            LA{localGeometry, localVariables, weights}.addJacobianAndResidualEntries(jacobian, residual);
        });
    }

private:
    template<typename LocalAssemblerFunc>
    void assemble_(const Variables& variables,
                   const LocalAssemblerFunc& assemble) const
    {
        auto fvGeometry = localView(*gridGeometry_);
        auto elemVariables = localView(variables);
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            fvGeometry.bind(element);
            elemVariables.bind(fvGeometry);
            assemble(fvGeometry, elemVariables);
        }
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;
};

} // namespace Dumux

#endif
