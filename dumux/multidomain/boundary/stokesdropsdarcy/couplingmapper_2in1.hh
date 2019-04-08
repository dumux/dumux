// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup StokesDropsDarcyCoupling
 * \copydoc Dumux::StokesDropsDarcyCouplingMapper
 */

#ifndef DUMUX_STOKES_DROPS_DARCY_COUPLINGMAPPER_HH
#define DUMUX_STOKES_DROPS_DARCY_COUPLINGMAPPER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <memory>
#include <unordered_map>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StokesDropsDarcyCoupling
 * \brief Coupling mapper for Stokes and Darcy domains with equal dimension87
 *
 *           and interface domain with lower dimension.
 */
template<class MDTraits>
class StokesDropsDarcyCouplingMapper
{
    using Scalar = typename MDTraits::Scalar;

public:
    static constexpr auto stokesCellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto stokesFaceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto cellCenterIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto faceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto stokesIdx = stokesCellCenterIdx;
    static constexpr auto interfaceIdx = typename MDTraits::template SubDomain<2>::Index();
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<3>::Index();

private:
    // obtain the type tags of the sub problems
    using StokesTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using InterfaceTypeTag = typename MDTraits::template SubDomain<2>::TypeTag;
    using DarcyTypeTag = typename MDTraits::template SubDomain<3>::TypeTag;

    struct ElementMapInfo // Stokes/Darcy to Interface
    {
        std::size_t eIdx; // interface element
        std::size_t scvfIdx; // own scvf
    };

    struct InterfaceElementMapInfo
    {
        std::size_t stokesEIdx;
        std::size_t stokesScvfIdx;
        std::size_t darcyEIdx;
        std::size_t darcyScvfIdx;
    };

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    using CouplingManager = GetPropType<StokesTypeTag, Properties::CouplingManager>;

    static_assert(GetPropType<SubDomainTypeTag<stokesIdx>, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::staggered,
                  "The free flow domain must use the staggered discretization");
    static_assert(GetPropType<SubDomainTypeTag<interfaceIdx>, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::cctpfa,
                  "The interface domain must use the CCTpfa discretization");
    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::cctpfa,
                  "The Darcy domain must use the CCTpfa discretization");

public:

    /*!
     * \brief Constructor
     */
    StokesDropsDarcyCouplingMapper(const CouplingManager& couplingManager) : couplingManager_(couplingManager) {}


    /*!
     * \brief
     */
    template<class Stencils>
    void computeCouplingMapsAndStencils(Stencils& stokesCellCenterToInterfaceStencils,
                                        Stencils& stokesFaceToInterfaceStencils,
                                        Stencils& darcyToInterfaceStencils,
                                        Stencils& interfaceToStokesCellCenterStencils,
                                        Stencils& interfaceToStokesFaceStencils,
                                        Stencils& interfaceToDarcyStencils)
    {
        const auto& stokesProblem = couplingManager_.problem(stokesIdx);
        const auto& interfaceProblem = couplingManager_.problem(interfaceIdx);
        const auto& darcyProblem = couplingManager_.problem(darcyIdx);

        const auto& stokesFvGridGeometry = stokesProblem.fvGridGeometry();
        const auto& interfaceFvGridGeometry = interfaceProblem.fvGridGeometry();
        const auto& darcyFvGridGeometry = darcyProblem.fvGridGeometry();

        isCoupledDarcyScvf_.resize(darcyFvGridGeometry.numScvf(), false);

        auto stokesFvGeometry = localView(stokesFvGridGeometry);
        auto interfaceFvGeometry = localView(interfaceFvGridGeometry);
        auto darcyFvGeometry = localView(darcyFvGridGeometry);

        const auto& stokesGridView = stokesFvGridGeometry.gridView();

        for(const auto& stokesElement : elements(stokesGridView))
        {
            stokesFvGeometry.bindElement(stokesElement);
            const auto stokesElementIdx = stokesFvGridGeometry.elementMapper().index(stokesElement);

            for(const auto& stokesScvf : scvfs(stokesFvGeometry))
            {
                // skip the DOF if it is not on the boundary
                if(!stokesScvf.boundary())
                    continue;

                // get element intersecting with the scvf center
                // for robustness, add epsilon in unit outer normal direction
                const auto eps = (stokesScvf.center() - stokesElement.geometry().center()).two_norm()*1e-8;
                auto globalPos = stokesScvf.center(); globalPos.axpy(eps, stokesScvf.unitOuterNormal());
                const auto darcyElementIdx = intersectingEntities(globalPos, darcyFvGridGeometry.boundingBoxTree());
                // TODO !!
//                const auto interfaceElementIdx = intersectingEntities(globalPos[0], interfaceFvGridGeometry.boundingBoxTree());
                const auto interfaceElementIdx = stokesElementIdx;

                // skip if no intersection was found
                if(darcyElementIdx.empty())
                    continue;

                // sanity check TODO
                if(darcyElementIdx.size() > 1)
                    DUNE_THROW(Dune::InvalidStateException, "Stokes face dof should only intersect with one interface element");

                const auto stokesCCDofIdx = stokesElementIdx;
                const auto stokesFaceDofIdx = stokesScvf.dofIndex();
                const auto darcyDofIdx = darcyElementIdx[0];
                const auto interfaceDofIdx = interfaceElementIdx;

                // fill stencils
                stokesCellCenterToInterfaceStencils[stokesCCDofIdx].push_back(interfaceDofIdx);
                stokesFaceToInterfaceStencils[stokesFaceDofIdx].push_back(interfaceDofIdx);
                darcyToInterfaceStencils[darcyDofIdx].push_back(interfaceDofIdx);
                interfaceToStokesCellCenterStencils[interfaceDofIdx].push_back(stokesCCDofIdx);
                interfaceToStokesFaceStencils[interfaceDofIdx].push_back(stokesFaceDofIdx);
                interfaceToDarcyStencils[interfaceDofIdx].push_back(darcyDofIdx);

                const auto& darcyElement = darcyFvGridGeometry.element(darcyElementIdx[0]);
                darcyFvGeometry.bindElement(darcyElement);

                // find the corresponding Darcy sub control volume face
                for(const auto& darcyScvf : scvfs(darcyFvGeometry))
                {
                    const Scalar distance = (darcyScvf.center() - stokesScvf.center()).two_norm();

                    if(distance < eps)
                    {
                        isCoupledDarcyScvf_[darcyScvf.index()] = true;

                        darcyElementToInterfaceElementMap_[darcyElementIdx[0]].push_back({interfaceElementIdx, darcyScvf.index()});
                        stokesElementToInterfaceElementMap_[stokesElementIdx].push_back({interfaceElementIdx, stokesScvf.index()});
                        interfaceElementMap_[interfaceElementIdx].push_back({stokesElementIdx, stokesScvf.index(), darcyElementIdx[0], darcyScvf.index()});

//                        // print maps and check all indices, dofs, ...
//                        std::cout << "** couplingmapper: interfaceElementMap for element " << interfaceElementIdx
//                                  << ": coupled stokes element = " << interfaceElementMap_[interfaceElementIdx][0].stokesEIdx
//                                  << ", coupled stokes scvf = " << interfaceElementMap_[interfaceElementIdx][0].stokesScvfIdx
//                                  << ", coupled darcy element = " << interfaceElementMap_[interfaceElementIdx][0].darcyEIdx
//                                  << ", coupled darcy scvf = " << interfaceElementMap_[interfaceElementIdx][0].darcyScvfIdx
//                                  << std::endl;

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
     * \brief A map that returns all Darcy elements coupled to a Stokes element
     */
    const auto& stokesElementToInterfaceElementMap() const
    {
        return stokesElementToInterfaceElementMap_;
    }

    /*!
     * \brief A map that returns information on all Darcy and Stokes elements coupled to an interface element
     */
    const auto& interfaceElementMap() const
    {
        return interfaceElementMap_;
    }

    /*!
     * \brief A map that returns information on all Darcy and Stokes elements coupled to an interface element
     */
    const auto& interfaceElementMap(const size_t interfaceElementIdx) const
    {
        return interfaceElementMap_.at(interfaceElementIdx);
    }

    /*!
     * \brief A map that returns all Stokes elements coupled to a Darcy element
     */
    const auto& darcyElementToInterfaceElementMap() const
    {
        return darcyElementToInterfaceElementMap_;
    }


private:
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> stokesElementToInterfaceElementMap_;
    std::unordered_map<std::size_t, std::vector<InterfaceElementMapInfo>> interfaceElementMap_;
    std::unordered_map<std::size_t, std::vector<ElementMapInfo>> darcyElementToInterfaceElementMap_;

    std::vector<bool> isCoupledDarcyScvf_;

    const CouplingManager& couplingManager_;
};

} // end namespace Dumux

#endif
