// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FoepplVonKarmanPlate
 * \brief A piece of the boundary whose shear gradient potential is one unknown constant
 *
 * Prescribing the shear gradient potential on a connected piece of the boundary drops the
 * compatibility rows of its nodes and with them the net flux of
 * \f$ \nabla w - \boldsymbol{\theta} \f$ through the piece. On one such piece the value is the
 * gauge; on every further one, such as a clamped hole in a clamped plate or a second
 * clamped edge joined to the first by free edges, it is an unknown constant.
 * The two classes here impose it: the compatibility rows of the piece's nodes are
 * replaced by \f$ \varphi_j = \varphi_0 \f$ and, at one master node, by the sum of all of
 * them. The coupling manager adds the entries the master row needs to the sparsity
 * pattern and to the coupling stencils, and the assembler performs the row operations
 * after every assembly.
 */
#ifndef DUMUX_FOEPPL_VON_KARMAN_PLATE_HOLE_GAUGE_HH
#define DUMUX_FOEPPL_VON_KARMAN_PLATE_HOLE_GAUGE_HH

#include <array>
#include <set>
#include <tuple>
#include <vector>
#include <cstddef>

#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/solidmechanics/plate/foeppl_von_karman/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Coupling manager of the Föppl-von Kármán plate with one piece of the boundary of
 *        the deformation domain carrying a single unknown value of the potential
 */
template<class MDTraits>
class HoleGaugeCouplingManager : public FoepplVonKarmanPlateCouplingManager<MDTraits>
{
    using ParentType = FoepplVonKarmanPlateCouplingManager<MDTraits>;
    static constexpr std::size_t numDomains = MDTraits::numSubDomains;
    static_assert(numDomains == 3, "The hole gauge is written for the three-domain plate");
    template<std::size_t id> using GridView = typename MDTraits::template SubDomain<id>::GridGeometry::GridView;
    template<std::size_t id> using Stencil = std::vector<typename IndexTraits<GridView<id>>::GridIndex>;

public:
    using ParentType::ParentType;
    static constexpr auto deformationIdx = ParentType::deformationIdx;

    /*!
     * \brief Declare the nodes of the piece, the first of which is the master
     *
     * The master row collects the compatibility rows of all nodes, so it couples to
     * every degree of freedom of every element that contains a node of the piece.
     * Call after init(), since the stencils of the sub-problems are needed, and before
     * the assembler is constructed, which builds the sparsity pattern.
     */
    template<class GridGeometry>
    void setHole(const GridGeometry& gridGeometry, const std::vector<std::size_t>& holeDofs)
    {
        hole_ = holeDofs;
        std::vector<bool> isHole(gridGeometry.numDofs(), false);
        for (const auto dof : hole_)
            isHole[dof] = true;

        const auto numElements = gridGeometry.gridView().size(0);
        isMasterElement_.assign(numElements, false);
        std::array<std::set<std::size_t>, numDomains> unions;
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bindElement(element);
            bool touchesHole = false, containsMaster = false;
            for (const auto& scv : scvs(fvGeometry))
            {
                touchesHole = touchesHole || isHole[scv.dofIndex()];
                containsMaster = containsMaster || scv.dofIndex() == hole_[0];
            }
            const auto eIdx = gridGeometry.elementMapper().index(element);
            isMasterElement_[eIdx] = containsMaster;
            if (!touchesHole)
                continue;
            for (const auto& scv : scvs(fvGeometry))
                unions[deformationIdx].insert(scv.dofIndex());
            Dune::Hybrid::forEach(std::make_index_sequence<numDomains>{}, [&](auto domainJ)
            {
                if constexpr (domainJ != deformationIdx)
                    for (const auto dof : this->ParentType::couplingStencil(Dune::index_constant<deformationIdx>{}, element, domainJ))
                        unions[domainJ].insert(dof);
            });
        }
        Dune::Hybrid::forEach(std::make_index_sequence<numDomains>{}, [&](auto j)
        {
            std::get<j>(unionStencils_).assign(unions[j].begin(), unions[j].end());
        });
    }

    const std::vector<std::size_t>& hole() const { return hole_; }

    template<std::size_t i, class Element, std::size_t j>
    const Stencil<j>& couplingStencil(Dune::index_constant<i> domainI,
                                      const Element& element,
                                      Dune::index_constant<j> domainJ) const
    {
        if constexpr (i == deformationIdx)
            if (!hole_.empty()
                && isMasterElement_[this->problem(domainI).gridGeometry().elementMapper().index(element)])
                return std::get<j>(unionStencils_);
        return ParentType::couplingStencil(domainI, element, domainJ);
    }

    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {
        if constexpr (id == deformationIdx)
        {
            if (hole_.empty())
                return;
            for (const auto dof : std::get<deformationIdx>(unionStencils_))
                pattern.add(hole_[0], dof);
            for (const auto dof : hole_)
                pattern.add(dof, hole_[0]);
        }
    }

private:
    std::vector<std::size_t> hole_;
    std::vector<bool> isMasterElement_;
    std::tuple<Stencil<0>, Stencil<1>, Stencil<2>> unionStencils_;
};

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief Assembler that replaces the compatibility rows of a piece of the boundary
 *        by the constraint of a single potential value and the summed row
 */
template<class BaseAssembler>
class HoleGaugeAssembler : public BaseAssembler
{
    static constexpr auto deformationIdx = BaseAssembler::CouplingManager::deformationIdx;

public:
    using SolutionVector = typename BaseAssembler::SolutionVector;
    using BaseAssembler::BaseAssembler;

    //! the nodes of the piece (master first), the potential's index and the row it replaces
    void setHole(std::vector<std::size_t> holeDofs, int potentialIdx, int compatibilityEqIdx)
    {
        hole_ = std::move(holeDofs);
        pvIdx_ = potentialIdx;
        eqIdx_ = compatibilityEqIdx;
    }

    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        BaseAssembler::assembleJacobianAndResidual(curSol);
        constrainJacobian_();
        constrainResidual_(curSol);
    }

    void assembleResidual(const SolutionVector& curSol)
    {
        BaseAssembler::assembleResidual(curSol);
        constrainResidual_(curSol);
    }

private:
    void constrainJacobian_()
    {
        if (hole_.empty())
            return;
        const auto master = hole_[0];
        auto& jacRow = this->jacobian()[Dune::index_constant<deformationIdx>{}];
        Dune::Hybrid::forEach(jacRow, [&](auto& block)
        {
            for (std::size_t n = 1; n < hole_.size(); ++n)
            {
                const auto j = hole_[n];
                for (auto colIt = block[j].begin(); colIt != block[j].end(); ++colIt)
                {
                    block[master][colIt.index()][eqIdx_] += (*colIt)[eqIdx_];
                    (*colIt)[eqIdx_] = 0.0;
                }
            }
        });
        auto& diagonal = jacRow[Dune::index_constant<deformationIdx>{}];
        for (std::size_t n = 1; n < hole_.size(); ++n)
        {
            diagonal[hole_[n]][hole_[n]][eqIdx_][pvIdx_] = 1.0;
            diagonal[hole_[n]][master][eqIdx_][pvIdx_] = -1.0;
        }
    }

    void constrainResidual_(const SolutionVector& curSol)
    {
        if (hole_.empty())
            return;
        const auto master = hole_[0];
        auto& residual = this->residual()[Dune::index_constant<deformationIdx>{}];
        const auto& sol = curSol[Dune::index_constant<deformationIdx>{}];
        for (std::size_t n = 1; n < hole_.size(); ++n)
        {
            residual[master][eqIdx_] += residual[hole_[n]][eqIdx_];
            residual[hole_[n]][eqIdx_] = sol[hole_[n]][pvIdx_] - sol[master][pvIdx_];
        }
    }

    std::vector<std::size_t> hole_;
    int pvIdx_ = 0, eqIdx_ = 1;
};

} // end namespace Dumux

#endif
