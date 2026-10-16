// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Constraints related to Dirichlet boundaries
 */
#ifndef DUMUX_DIRICHLET_CONSTRAINTS_HH
#define DUMUX_DIRICHLET_CONSTRAINTS_HH

#include <array>
#include <bitset>
#include <concepts>
#include <cstddef>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/multimapperview.hh>

namespace Dumux {

template<class Info, class Values>
struct ConstraintData
{
    using ConstraintInfo = Info;
    using ConstraintValues = Values;

    const ConstraintInfo& constraintInfo() const
    { return info_; }

    const ConstraintValues& values() const
    { return values_; }

    ConstraintInfo info_;
    ConstraintValues values_;
};

template<class DirichletConstraintInfo, class DirichletValues, class IndexType>
struct DirichletConstraintData : public ConstraintData<DirichletConstraintInfo, DirichletValues>
{
    using GridIndexType = IndexType;

    GridIndexType dofIndex() const
    { return dofIdx_; }

    GridIndexType dofIdx_;
};

namespace Detail {

//! The Dirichlet constraint of one equation of a dof: whether it is set, the constrained primary variable and its value
struct ConstrainedEquation
{
    bool isConstrained = false;
    int pvIdx = 0;
    double value = 0.0;
};

/*!
 * \brief Sends the constrained equations of the dofs on process-border entities
 *        and collects the received ones that are not constrained locally
 */
template<class GridView, class DofMapper, int numEq>
class ConstraintExchange
: public Dune::CommDataHandleIF<ConstraintExchange<GridView, DofMapper, numEq>, double>
{
    static constexpr int dim = GridView::dimension;
    using Equations = std::array<ConstrainedEquation, numEq>;
    // per dof and equation: whether it is constrained, its primary variable and value
    static constexpr std::size_t itemsPerDof = 3*numEq;

public:
    ConstraintExchange(const DofMapper& dofMapper,
                       const std::unordered_map<std::size_t, Equations>& local,
                       std::map<std::size_t, Equations>& received)
    : dofMapper_(dofMapper)
    , local_(local)
    , received_(received)
    {
        using Grid = typename GridView::Grid;
        Dune::Hybrid::forEach(std::make_index_sequence<dim+1>{}, [&](auto codim)
        {
            if constexpr (Dune::Capabilities::canCommunicate<Grid, codim>::v)
                communicates_[codim] = !dofMapper_.types(codim).empty();
        });
    }

    bool contains(int, int codim) const
    { return codim >= 0 && codim <= dim && communicates_[codim]; }

    bool fixedSize(int, int) const
    { return false; }

    template<class Entity>
    std::size_t size(const Entity& entity) const
    { return asMultiMapper(dofMapper_).indices(entity).size()*itemsPerDof; }

    template<class MessageBuffer, class Entity>
    void gather(MessageBuffer& buffer, const Entity& entity) const
    {
        for (const auto dofIdx : asMultiMapper(dofMapper_).indices(entity))
        {
            const auto it = local_.find(dofIdx);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                const auto equation = it != local_.end() ? it->second[eqIdx] : ConstrainedEquation{};
                buffer.write(equation.isConstrained ? 1.0 : 0.0);
                buffer.write(static_cast<double>(equation.pvIdx));
                buffer.write(equation.value);
            }
        }
    }

    template<class MessageBuffer, class Entity>
    void scatter(MessageBuffer& buffer, const Entity& entity, std::size_t)
    {
        for (const auto dofIdx : asMultiMapper(dofMapper_).indices(entity))
        {
            const auto it = local_.find(dofIdx);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                double isConstrained, pvIdx, value;
                buffer.read(isConstrained);
                buffer.read(pvIdx);
                buffer.read(value);

                const bool isLocal = it != local_.end() && it->second[eqIdx].isConstrained;
                if (isConstrained > 0.5 && !isLocal && !received_[dofIdx][eqIdx].isConstrained)
                    received_[dofIdx][eqIdx] = ConstrainedEquation{true, static_cast<int>(pvIdx), value};
            }
        }
    }

private:
    const DofMapper& dofMapper_;
    const std::unordered_map<std::size_t, Equations>& local_;
    std::map<std::size_t, Equations>& received_;
    std::bitset<dim+1> communicates_;
};

} // end namespace Detail

/*!
 * \ingroup Discretization
 * \brief Make a list of Dirichlet constraints consistent between the processes of a parallel run
 *
 * A process knows the boundary faces of its own (and ghost) elements. A dof on the domain boundary can belong to
 * elements of a process none of which has a boundary face containing the dof, so constraints collected from boundary
 * faces miss that dof on that process, and its rows assembled there are not constrained. For each dof on a border
 * entity (shared by the interior partitions of several processes), this function adds a constraint of the equations
 * that other processes constrain and this process does not. Dofs on overlap entities are not exchanged: with an overlap
 * of at least one element layer, a process has all elements around its interior and border dofs. The constraints of
 * a dof may be split over several entries, e.g. one per boundary face. An equation constrained on several processes
 * keeps the local constraint (the values are expected to agree). The list does not change in a sequential run.
 *
 * \param gridDiscretization the grid discretization whose dof mapper numbers the constrained dofs
 * \param constraints the constraints of this process (DirichletConstraintData), extended in place
 */
template<class GridDiscretization, class Constraint>
void makeConstraintsConsistent(const GridDiscretization& gridDiscretization, std::vector<Constraint>& constraints)
{
    const auto& gridView = gridDiscretization.gridView();
    if (gridView.comm().size() <= 1)
        return;

    constexpr int numEq = Constraint::ConstraintValues::dimension;
    using Equations = std::array<Detail::ConstrainedEquation, numEq>;

    // a later entry constraining the same equation overrides an earlier one, as in the assembly
    std::unordered_map<std::size_t, Equations> local;
    for (const auto& constraint : constraints)
    {
        const auto& info = constraint.constraintInfo();
        auto& equations = local[constraint.dofIndex()];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (info.isConstraintEquation(eqIdx))
            {
                const int pvIdx = info.eqToPriVarIndex(eqIdx);
                equations[eqIdx] = Detail::ConstrainedEquation{true, pvIdx, static_cast<double>(constraint.values()[pvIdx])};
            }
        }
    }

    const auto& dofMapper = gridDiscretization.dofMapper();
    using GridView = std::decay_t<decltype(gridView)>;
    using DofMapper = std::decay_t<decltype(dofMapper)>;
    std::map<std::size_t, Equations> received;
    Detail::ConstraintExchange<GridView, DofMapper, numEq> exchange(dofMapper, local, received);
    gridView.communicate(exchange, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

    for (const auto& [dofIdx, equations] : received)
    {
        typename Constraint::ConstraintInfo info;
        typename Constraint::ConstraintValues values(0.0);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (equations[eqIdx].isConstrained)
            {
                info.set(equations[eqIdx].pvIdx, eqIdx);
                values[equations[eqIdx].pvIdx] = equations[eqIdx].value;
            }
        }

        constraints.push_back(Constraint{std::move(info), std::move(values), static_cast<typename Constraint::GridIndexType>(dofIdx)});
    }
}

} // end namespace Dumux

#endif
