// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPoreNetworkCoupling
 * \copydoc Dumux::StaggeredFreeFlowPoreNetworkCouplingMapper
 */

#ifndef DUMUX_MULTIDOMAIN_FREEFLOW_POROUSMEDIUM_COUPLINGMAPPER_HH
#define DUMUX_MULTIDOMAIN_FREEFLOW_POROUSMEDIUM_COUPLINGMAPPER_HH

#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {

/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief Coupling mapper for staggered free-flow and pore-network models.
 */
// template<class MDTraits, class CouplingManager>
class StaggeredFreeFlowPoreNetworkCouplingMapper
{
    using MapType = std::unordered_map<std::size_t, std::vector<std::size_t>>;

public:
    /*!
     * \brief Main update routine
     */
    template<class FreeFlowMomentumGridGeometry, class FreeFlowMassGridGeometry, class PoreNetworkGridGeometry>
    void update(const FreeFlowMomentumGridGeometry& ffMomentumGridGeometry,
                const FreeFlowMassGridGeometry& ffMassGridGeometry,
                const PoreNetworkGridGeometry& pnmGridGeometry)
    {
        clear_();
        resize_(ffMomentumGridGeometry, pnmGridGeometry);
        Dune::Timer watch;
        std::cout << "Initializing the coupling map..." << std::endl;

        auto ffFvGeometry = localView(ffMomentumGridGeometry);
        auto pnmFvGeometry = localView(pnmGridGeometry);

        using GlobalPosition = typename FreeFlowMomentumGridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

        for (const auto& pnmElement : elements(pnmGridGeometry.gridView()))
        {
            const auto pnmElementIdx = pnmGridGeometry.elementMapper().index(pnmElement);
            pnmFvGeometry.bindElement(pnmElement);
            for (const auto& pnmScv : scvs(pnmFvGeometry))
            {
                // skip the dof if it is not on the boundary
                if (!pnmGridGeometry.dofOnBoundary(pnmScv.dofIndex()))
                    continue;

                // get the intersection bulk element
                const auto pnmPos = pnmScv.dofPosition();
                const auto pnmDofIdx = pnmScv.dofIndex();
                const auto& otherPNMScv = pnmFvGeometry.scv(1 - pnmScv.indexInElement());
                const auto otherPNMScvDofIdx = otherPNMScv.dofIndex();

                // check for intersections, skip if no intersection was found
                const auto directlyCoupledFreeFlowElements = intersectingEntities(pnmPos, ffMomentumGridGeometry.boundingBoxTree());
                if (directlyCoupledFreeFlowElements.empty())
                    continue;
                else
                    isCoupledPNMDof_[pnmDofIdx] = true;

                // determine the normal direction of the local coupling interface heuristically:
                // find all element intersections touching the pore and take the normal which
                // occurs most frequently
                const std::size_t couplingNormalDirectionIndex = [&]
                {
                    using Key = std::pair<std::size_t, bool>;
                    std::map<Key, std::size_t> result;
                    for (const auto eIdx : directlyCoupledFreeFlowElements)
                    {
                        for (const auto& intersection : intersections(ffMomentumGridGeometry.gridView(), ffMomentumGridGeometry.element(eIdx)))
                        {
                            if (intersectsPointGeometry(pnmPos, intersection.geometry()))
                            {
                                const auto& normal = intersection.centerUnitOuterNormal();
                                const auto normalAxis = Dumux::normalAxis(normal);
                                const Key key = std::make_pair(normalAxis, std::signbit(normal[normalAxis]));
                                ++result[key];
                            }
                        }
                    }

                    // TODO how to properly handle this corner (literally) case
                    if (directlyCoupledFreeFlowElements.size() == 1 && result.size() > 1)
                        DUNE_THROW(Dune::InvalidStateException, "Pore may not intersect with faces of different orientation when coupled to only one element");

                    return std::max_element(result.begin(), result.end(), [](const auto& x, const auto& y) { return x.second < y.second;})->first.first;
                }();

                using Scalar = typename FreeFlowMomentumGridGeometry::GridView::ctype;
                const Scalar couplingPoreRadius = pnmGridGeometry.poreInscribedRadius(pnmDofIdx);
                GlobalPosition lowerLeft = pnmPos - GlobalPosition(couplingPoreRadius);
                lowerLeft[couplingNormalDirectionIndex] = pnmPos[couplingNormalDirectionIndex];
                GlobalPosition upperRight = pnmPos + GlobalPosition(couplingPoreRadius);
                upperRight[couplingNormalDirectionIndex] = pnmPos[couplingNormalDirectionIndex];

                auto axes = std::move(std::bitset<FreeFlowMomentumGridGeometry::Grid::dimensionworld>{}.set());
                axes.set(couplingNormalDirectionIndex, false);

                using PoreIntersectionGeometryType = Dune::AxisAlignedCubeGeometry<Scalar,
                                                                                   FreeFlowMomentumGridGeometry::GridView::dimension-1,
                                                                                   FreeFlowMomentumGridGeometry::GridView::dimensionworld>;

                PoreIntersectionGeometryType poreIntersectionGeometry(lowerLeft, upperRight, axes);
                const auto allCoupledFreeFlowElements = intersectingEntities(std::move(poreIntersectionGeometry), ffMomentumGridGeometry.boundingBoxTree());


                for (const auto& ffElementInfo : allCoupledFreeFlowElements)
                {
                    const auto freeFlowElementIndex = ffElementInfo.second();
                    pnmElementToFreeFlowElementsMap_[pnmElementIdx].push_back(freeFlowElementIndex);
                    freeFlowElementToPNMElementMap_[freeFlowElementIndex] = pnmElementIdx;

                    pnmToFreeFlowMassStencils_[pnmElementIdx].push_back(freeFlowElementIndex);
                    freeFlowMassToPNMStencils_[freeFlowElementIndex].push_back(pnmDofIdx);

                    ffFvGeometry.bindElement(ffMomentumGridGeometry.element(freeFlowElementIndex));
                    const auto coupledFreeFlowMomentumDofIndices = coupledFFMomentumDofs_(ffFvGeometry, ffMassGridGeometry, pnmPos, couplingPoreRadius, couplingNormalDirectionIndex);

                    pnmToFreeFlowMomentumStencils_[pnmElementIdx].push_back(coupledFreeFlowMomentumDofIndices.coupledFrontalDof);
                    freeFlowMomentumToPNMStencils_[coupledFreeFlowMomentumDofIndices.coupledFrontalDof].push_back(pnmDofIdx);
                    freeFlowMomentumToPNMStencils_[coupledFreeFlowMomentumDofIndices.coupledFrontalDof].push_back(otherPNMScvDofIdx);

                    isCoupledFreeFlowMomentumDof_[coupledFreeFlowMomentumDofIndices.coupledFrontalDof] = true;
                    isCoupledFreeFlowMomentumDofOnInterface_[coupledFreeFlowMomentumDofIndices.coupledFrontalDof] = true;

                    // treat the coupled ff dofs not directly on interface
                    for (const auto ffDofIdx : coupledFreeFlowMomentumDofIndices.coupledLateralDofs)
                    {
                        freeFlowMomentumToPNMStencils_[ffDofIdx].push_back(pnmDofIdx);
                        freeFlowMomentumToPNMStencils_[ffDofIdx].push_back(otherPNMScvDofIdx);
                        isCoupledFreeFlowMomentumDof_[ffDofIdx] = true;
                    }
                }
            }
        }

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     * \param eIdxI the index of the coupled element of domain í
     */
    const std::vector<std::size_t>& poreNetworkToFreeFlowMomentumCouplingStencil(const std::size_t eIdxI) const
    {
        if (isCoupledPoreNetworkElement(eIdxI))
            return pnmToFreeFlowMomentumStencils_.at(eIdxI);
        else
            return emptyStencil_;
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     * \param eIdxI the index of the coupled element of domain í
     */
    const std::vector<std::size_t>& poreNetworkToFreeFlowMassCouplingStencil(const std::size_t eIdxI) const
    {
        if (isCoupledPoreNetworkElement(eIdxI))
            return pnmToFreeFlowMassStencils_.at(eIdxI);
        else
            return emptyStencil_;
    }

        /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     * \param eIdxI the index of the coupled element of domain i
     */
    const std::vector<std::size_t>& freeFlowMassToPoreNetworkCouplingStencil(const std::size_t eIdxI) const
    {
        if (isCoupledFreeFlowElement(eIdxI))
            return freeFlowMassToPNMStencils_.at(eIdxI);
        else
            return emptyStencil_;
    }

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     * \param dofIndex the degree of freedom index
     */
    const std::vector<std::size_t>& freeFlowMomentumToPoreNetworkCouplingStencil(const std::size_t dofIndex) const
    {
        if (isCoupledFreeFlowMomentumDof(dofIndex))
            return freeFlowMomentumToPNMStencils_.at(dofIndex);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Return if an element residual with index eIdx of domain i is coupled to domain j
     */
    bool isCoupledFreeFlowElement(std::size_t eIdx) const
    {
        return static_cast<bool>(freeFlowElementToPNMElementMap_.count(eIdx));
    }

    /*!
     * \brief Return if an element residual with index eIdx of domain i is coupled to domain j
     */
    bool isCoupledFreeFlowMomentumDof(std::size_t dofIdx) const
    {
        return isCoupledFreeFlowMomentumDof_[dofIdx];
    }

    /*!
     * \brief Return if an element residual with index eIdx of domain i is coupled to domain j
     */
    bool isCoupledPoreNetworkElement(std::size_t eIdx) const
    {
        return static_cast<bool>(pnmElementToFreeFlowElementsMap_.count(eIdx));
    }

    /*!
     * \brief Return if an element residual with index eIdx of domain i is coupled to domain j
     */
    bool isCoupledPoreNetworkDof(std::size_t dofIdx) const
    {
        return isCoupledPNMDof_[dofIdx];
    }

    bool isCoupledFreeFlowMomentumScvf(std::size_t scvfIdx) const
    {
        return isCoupledFrontalFreeFlowMomentumScvf_.count(scvfIdx);
    }

    bool isCoupledFreeFlowMomentumLateralScvf(std::size_t scvfIdx) const
    {
        return isCoupledLateralFreeFlowMomentumScvf_.count(scvfIdx);
    }

    bool isCoupledFreeFlowMassScvf(std::size_t scvfIdx) const
    {
        return isCoupledFreeFlowMassScvf_.count(scvfIdx);
    }

    const auto& pnmElementToFreeFlowElementsMap() const
    { return pnmElementToFreeFlowElementsMap_;}

    const auto& freeFlowElementToPNMElementMap() const
    { return freeFlowElementToPNMElementMap_; }

private:

    void clear_()
    {
        pnmElementToFreeFlowElementsMap_.clear();
        freeFlowElementToPNMElementMap_.clear();
        isCoupledPNMDof_.clear();
        isCoupledFreeFlowMomentumDof_.clear();
        isCoupledFreeFlowMomentumDofOnInterface_.clear();
        pnmToFreeFlowMassStencils_.clear();
        pnmToFreeFlowMomentumStencils_.clear();
        freeFlowMassToPNMStencils_.clear();
        freeFlowMomentumToPNMStencils_.clear();
    }

    template<class FreeFlowMomentumGridGeometry, class PoreNetworkGridGeometry>
    void resize_(const FreeFlowMomentumGridGeometry& ffMomentumGridGeometry,
                 const PoreNetworkGridGeometry& pnmGridGeometry)
    {
        const auto numPNMDofs = pnmGridGeometry.numDofs();
        const auto numFreeFlowMomentumDofs = ffMomentumGridGeometry.numDofs();
        isCoupledPNMDof_.resize(numPNMDofs, false);
        isCoupledFreeFlowMomentumDof_.resize(numFreeFlowMomentumDofs, false);
        isCoupledFreeFlowMomentumDofOnInterface_.resize(numFreeFlowMomentumDofs, false);
    }

     //! get the free-flow momentum dofs that are coupled to the PNM dofs
    template<class FVElementGeometry, class FreeFlowMassGridGeometry, class GlobalPosition, class Scalar>
    auto coupledFFMomentumDofs_(const FVElementGeometry& fvGeometry,
                                const FreeFlowMassGridGeometry& ffMassGridGeometry,
                                const GlobalPosition& pnmPos,
                                const Scalar couplingPoreRadius,
                                const int couplingInterfaceDirectionIdx)
    {

        struct Result
        {
            Dune::ReservedVector<std::size_t, FVElementGeometry::maxNumElementScvs> coupledLateralDofs;
            std::size_t coupledFrontalDof;
        } result;


        using std::abs;
        for (const auto& scv : scvs(fvGeometry))
        {
            const Scalar eps = diameter(fvGeometry.geometry(scv))*1e-6; // TODO
            assert(eps < couplingPoreRadius);

            if (scv.dofAxis() == couplingInterfaceDirectionIdx) // the free flow dofs that lie within the coupling interface
            {
                if (abs(scv.dofPosition()[couplingInterfaceDirectionIdx] - pnmPos[couplingInterfaceDirectionIdx]) < eps)
                {
                    result.coupledFrontalDof = scv.dofIndex();

                    // treat scvfs
                    for (const auto& scvf : scvfs(fvGeometry, scv))
                    {
                        // add lateral faces "standing" on coupling interface
                        if (scvf.isLateral() && !scvf.boundary())
                            isCoupledLateralFreeFlowMomentumScvf_[scvf.index()] = true;
                        else if (scvf.isFrontal() && scvf.boundary()) // add face lying on interface
                        {
                            isCoupledFrontalFreeFlowMomentumScvf_[scvf.index()] = true;

                            const auto& element = ffMassGridGeometry.element(fvGeometry.elementIndex()); // this local variable is needed to prevent a memory error
                            const auto ffMassFVGeometry = localView(ffMassGridGeometry).bindElement(element);
                            for (const auto& ffMassScvf : scvfs(ffMassFVGeometry))
                            {
                                if (abs(ffMassScvf.center()[couplingInterfaceDirectionIdx] - pnmPos[couplingInterfaceDirectionIdx]) < eps)
                                    isCoupledFreeFlowMassScvf_[ffMassScvf.index()] = true;
                            }
                        }
                    }
                }
            }
            else // the free flow dofs perpendicular to the coupling interface
            {
                bool isCoupledDof = false;

                for (int dimIdx = 0; dimIdx < GlobalPosition::dimension; ++dimIdx)
                {
                    if (dimIdx == couplingInterfaceDirectionIdx || scv.boundary())
                        continue;

                    isCoupledDof = abs(scv.dofPosition()[dimIdx] - pnmPos[dimIdx]) < couplingPoreRadius + eps;
                }

                if (isCoupledDof)
                {
                    result.coupledLateralDofs.push_back(scv.dofIndex());

                    // treat scvfs
                    for (const auto& scvf : scvfs(fvGeometry, scv))
                    {
                        if (scvf.isLateral() && scvf.boundary())
                        {
                            // add lateral scvfs lying on interface
                            if (abs(scvf.ipGlobal()[couplingInterfaceDirectionIdx] - pnmPos[couplingInterfaceDirectionIdx]) < eps)
                                isCoupledLateralFreeFlowMomentumScvf_[scvf.index()] = true;
                        }
                    }
                }
            }
        }

        return result;
    }

    std::vector<std::size_t> emptyStencil_;

    std::unordered_map<std::size_t, std::vector<std::size_t>> pnmElementToFreeFlowElementsMap_;
    std::unordered_map<std::size_t, std::size_t> freeFlowElementToPNMElementMap_;

    std::vector<bool> isCoupledPNMDof_;
    std::vector<bool> isCoupledFreeFlowMomentumDof_;
    std::vector<bool> isCoupledFreeFlowMomentumDofOnInterface_;

    std::unordered_map<std::size_t, bool> isCoupledLateralFreeFlowMomentumScvf_;
    std::unordered_map<std::size_t, bool> isCoupledFrontalFreeFlowMomentumScvf_;
    std::unordered_map<std::size_t, bool> isCoupledFreeFlowMassScvf_;

    MapType pnmToFreeFlowMassStencils_;
    MapType pnmToFreeFlowMomentumStencils_;
    MapType freeFlowMassToPNMStencils_;
    MapType freeFlowMomentumToPNMStencils_;
};

} // end namespace Dumux

#endif
