// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingMapper
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGMAPPER_HH
#define DUMUX_STOKES_DARCY_COUPLINGMAPPER_HH

#include <type_traits>
#include <unordered_map>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension.
 */
class StokesDarcyCouplingMapper
{
    struct ElementMapInfo
    {
        std::size_t eIdx;
        std::size_t scvfIdx;
        std::size_t flipScvfIdx;
    };

public:

    /*!
     * \brief Main update routine
     */
    template<class CouplingManager, class Stencils>
    void computeCouplingMapsAndStencils(const CouplingManager& couplingManager,
                                        Stencils& darcyToStokesCellCenterStencils,
                                        Stencils& darcyToStokesFaceStencils,
                                        Stencils& stokesCellCenterToDarcyStencils,
                                        Stencils& stokesFaceToDarcyStencils)
    {
        const auto& stokesFvGridGeometry = couplingManager.problem(CouplingManager::stokesIdx).gridGeometry();
        const auto& darcyFvGridGeometry = couplingManager.problem(CouplingManager::darcyIdx).gridGeometry();

        static_assert(std::decay_t<decltype(stokesFvGridGeometry)>::discMethod == DiscretizationMethods::staggered,
                      "The free flow domain must use the staggered discretization");

        static_assert(std::decay_t<decltype(darcyFvGridGeometry)>::discMethod == DiscretizationMethods::cctpfa,
                      "The Darcy domain must use the CCTpfa discretization");

        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.numScvf(), false);

        const auto& stokesGridView = stokesFvGridGeometry.gridView();
        auto stokesFvGeometry = localView(stokesFvGridGeometry);
        for (const auto& stokesElement : elements(stokesGridView))
        {
            stokesFvGeometry.bindElement(stokesElement);

            for (const auto& scvf : scvfs(stokesFvGeometry))
            {
                // skip the DOF if it is not on the boundary
                if  (!scvf.boundary())
                    continue;

                // get element intersecting with the scvf center
                // for robustness, add epsilon in unit outer normal direction
                const auto eps = (scvf.center() - stokesElement.geometry().center()).two_norm()*1e-8;
                auto globalPos = scvf.center(); globalPos.axpy(eps, scvf.unitOuterNormal());
                const auto darcyElementIdx = intersectingEntities(globalPos, darcyFvGridGeometry.boundingBoxTree());

                // skip if no intersection was found
                if (darcyElementIdx.empty())
                    continue;

                // sanity check
                if (darcyElementIdx.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one Darcy element");

                const auto stokesElementIdx = stokesFvGridGeometry.elementMapper().index(stokesElement);

                const auto darcyDofIdx = darcyElementIdx[0];

                stokesFaceToDarcyStencils[scvf.dofIndex()].push_back(darcyDofIdx);
                stokesCellCenterToDarcyStencils[stokesElementIdx].push_back(darcyDofIdx);

                darcyToStokesFaceStencils[darcyElementIdx[0]].push_back(scvf.dofIndex());
                darcyToStokesCellCenterStencils[darcyElementIdx[0]].push_back(stokesElementIdx);

                const auto& darcyElement = darcyFvGridGeometry.element(darcyElementIdx[0]);
                const auto darcyFvGeometry = localView(darcyFvGridGeometry).bindElement(darcyElement);

                // find the corresponding Darcy sub control volume face
                for (const auto& darcyScvf : scvfs(darcyFvGeometry))
                {
                    const auto distance = (darcyScvf.center() - scvf.center()).two_norm();

                    if (distance < eps)
                    {
                        isCoupledDarcyScvf_[darcyScvf.index()] = true;
                        darcyElementToStokesElementMap_[darcyElementIdx[0]].push_back({stokesElementIdx, scvf.index(), darcyScvf.index()});
                        stokesElementToDarcyElementMap_[stokesElementIdx].push_back({darcyElementIdx[0], darcyScvf.index(), scvf.index()});
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns whether a Darcy scvf is coupled to the other domain
     */
    bool isCoupledDarcyScvf(std::size_t darcyScvfIdx) const
    {
        return isCoupledDarcyScvf_[darcyScvfIdx];
    }

    /*!
     * \brief A map that returns all Stokes elements coupled to a Darcy element
     */
    const auto& darcyElementToStokesElementMap() const
    {
        return darcyElementToStokesElementMap_;
    }

    /*!
     * \brief A map that returns all Darcy elements coupled to a Stokes element
     */
    const auto& stokesElementToDarcyElementMap() const
    {
        return stokesElementToDarcyElementMap_;
    }

private:
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> darcyElementToStokesElementMap_;
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> stokesElementToDarcyElementMap_;

    std::vector<bool> isCoupledDarcyScvf_;
};

} // end namespace Dumux

#endif
