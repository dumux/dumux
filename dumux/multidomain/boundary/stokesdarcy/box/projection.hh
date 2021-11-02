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
 * \ingroup StokesDarcyCoupling
 */

#ifndef DUMUX_STOKES_DARCY_PROJECTION_HH
#define DUMUX_STOKES_DARCY_PROJECTION_HH

#include <type_traits>

#include <dumux/common/properties.hh>
#include "couplingmanager.hh"

namespace Dumux {

namespace Detail {
    enum class ProjectionMethod
    {
        L2Projection, AreaWeightedDofEvaluation
    };

    template <typename T>
    using ProjectionMethodDetector = decltype(std::declval<T>().projectionMethod());

    template<class T>
    static constexpr bool hasProjectionMethod()
    { return Dune::Std::is_detected<ProjectionMethodDetector, T>::value; }


    template<class T>
    static constexpr ProjectionMethod projectionMethod()
    {
        if constexpr (hasProjectionMethod<T>())
            return T::projectionMethod();
        else
            return ProjectionMethod::L2Projection;
    }
}

template<class MDTraits, class CouplingManager>
class Projection
{
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto freeFlowIdx = CouplingManager::freeFlowIdx;
    static constexpr auto porousMediumIdx = CouplingManager::porousMediumIdx;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    using SolutionVector = typename MDTraits::SolutionVector;

    using ProjectionMethod = Detail::ProjectionMethod;
    static constexpr auto projectionMethod = Detail::projectionMethod<MDTraits>();

public:

    Projection(const SolutionVector& sol) : sol_(sol)
    { }

    // calculate projection of pm solution needed for fluxes of the free-flow residual
    template<class Function>
    Scalar calculateProjection(const CouplingManager& couplingManager,
                               const SubControlVolumeFace<freeFlowIdx>& stokesScvf,
                               const Element<porousMediumIdx>& darcyElement,
                               const ElementVolumeVariables<porousMediumIdx>& darcyElemVolVars,
                               Function evalPriVar) const
    {
        Scalar projection = 0.0;
        auto domainI = Dune::index_constant<freeFlowIdx>();
        auto fvGeometry = localView(couplingManager.problem(porousMediumIdx).gridGeometry());
        auto elemVolVars = localView(darcyElemVolVars.gridVolVars());

        // integrate darcy pressure over each coupling facet and average
        for(const auto& couplingFacet : couplingFacets(domainI, couplingManager.couplingMapper(), stokesScvf.insideScvIdx(), stokesScvf.localFaceIdx()))
        {
            const auto darcyEIdxI = couplingManager.problem(porousMediumIdx).gridGeometry().elementMapper().index(darcyElement);
            const auto darcyEIdxJ = couplingFacet.pmEIdx;

            const auto& element = couplingManager.problem(porousMediumIdx).gridGeometry().boundingBoxTree().entitySet().entity(couplingFacet.pmEIdx);
            fvGeometry.bind(element);

            if(darcyEIdxI == darcyEIdxJ)
            {
                projection += calculateFacetIntegral(element, fvGeometry, fvGeometry.scvf(couplingFacet.pmScvfIdx), darcyElemVolVars, couplingFacet.geometry, evalPriVar);
            }
            else
            {
                elemVolVars.bind(element, fvGeometry, sol_[porousMediumIdx]);
                projection += calculateFacetIntegral(element, fvGeometry, fvGeometry.scvf(couplingFacet.pmScvfIdx), elemVolVars, couplingFacet.geometry, evalPriVar);
            }

        }

        projection /= stokesScvf.area();

        return projection;
    }

    // calculate projection of pm solution needed for fluxes of ff residual
    template<class Function>
    Scalar calculateProjection(const CouplingManager& couplingManager,
                               const Element<freeFlowIdx>& element,
                               const SubControlVolumeFace<freeFlowIdx>& scvf,
                               Function evalPriVar) const
    {
        Scalar projection = 0.0;

        // integrate darcy pressure over each coupling facet and average
        for (const auto& data : couplingManager.stokesCouplingContext())
        {
            //ToDo Is this if really necessary?
            if (scvf.index() == data.stokesScvfIdx)
            {
                const auto& elemVolVars = data.elementVolVars;
                const auto& darcyScvf = data.fvGeometry.scvf(data.darcyScvfIdx);
                const auto& couplingFacet = couplingManager.couplingMapper().couplingFacet(data.facetIdx);
                projection += calculateFacetIntegral(data.element, data.fvGeometry, darcyScvf, elemVolVars, couplingFacet.geometry, evalPriVar);
            }
        }

        projection /= scvf.area();

        return projection;
    }

    template<class CouplingFacetGeometry, class ElementVolumeVariables, class Function>
    static Scalar calculateFacetIntegral(const Element<porousMediumIdx>& element,
                                         const FVElementGeometry<porousMediumIdx>& fvGeometry,
                                         const SubControlVolumeFace<porousMediumIdx>& scvf,
                                         const ElementVolumeVariables& elemVolVars,
                                         const CouplingFacetGeometry& facetGeometry,
                                         Function evalPriVar)
    {
        Scalar facetProjection = 0.0;
        if constexpr (projectionMethod == ProjectionMethod::L2Projection)
        {
            const auto& localBasis = fvGeometry.feLocalBasis();

            // do second order integration as box provides linear functions
            static constexpr int darcyDim = GridGeometry<porousMediumIdx>::GridView::dimension;
            const auto& rule = Dune::QuadratureRules<Scalar, darcyDim-1>::rule(facetGeometry.type(), 2);
            for (const auto& qp : rule)
            {
                const auto& ipLocal = qp.position();
                const auto& ipGlobal = facetGeometry.global(ipLocal);
                const auto& ipElementLocal = element.geometry().local(ipGlobal);

                std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
                localBasis.evaluateFunction(ipElementLocal, shapeValues);

                Scalar value = 0.0;
                for (const auto& scv : scvs(fvGeometry))
                    value += evalPriVar(elemVolVars[scv.localDofIndex()])*shapeValues[scv.indexInElement()][0];

                facetProjection += value*facetGeometry.integrationElement(qp.position())*qp.weight();
            }
        }
        else if constexpr (projectionMethod == ProjectionMethod::AreaWeightedDofEvaluation)
        {
            facetProjection = facetGeometry.volume()*evalPriVar(elemVolVars[fvGeometry.scv(scvf.insideScvIdx()).localDofIndex()]);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Unkown projection method!");
        }

        return facetProjection;
    }

private:

const SolutionVector& sol_;

};

} // end namespace Dumux

#endif // DUMUX_STOKES_DARCY_PROJECTION_HH
