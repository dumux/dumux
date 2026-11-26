// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Tissue and network problems of the embedded coupling test on distributed grids
 *
 * The exchange is a membrane flux 2 pi R c (p_network - p_tissue) with the radius R of the network
 * element and a coefficient c of the network element given to the coupling manager as element data. It
 * is scaled with (1 + d) by the average distance d to the bulk element where the coupling mode provides
 * it, and with the flux scaling factor for the kernel coupling. The coupling terms can be switched off to
 * separate them from the rest of the residual.
 */
#ifndef DUMUX_TEST_EMBEDDED_PARALLEL_PROBLEMS_HH
#define DUMUX_TEST_EMBEDDED_PARALLEL_PROBLEMS_HH

#include <cmath>
#include <vector>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

namespace Dumux {

//! The exchange per unit length of the network for the point source id
template<class CouplingManager>
auto embeddedParallelExchange(const CouplingManager& couplingManager, std::size_t id)
{
    const auto p3D = couplingManager.bulkPriVars(id)[0];
    const auto p1D = couplingManager.lowDimPriVars(id)[0];
    auto exchange = 2.0*M_PI*couplingManager.radius(id)*couplingManager.lowDimElementData(id, 0)*(p1D - p3D);

    if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
        exchange *= couplingManager.fluxScalingFactor(id);
    if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel
                  || CouplingManager::couplingMode == Embedded1d3dCouplingMode::line
                  || CouplingManager::couplingMode == Embedded1d3dCouplingMode::average)
        exchange *= 1.0 + couplingManager.averageDistance(id);

    return exchange;
}

template<class TypeTag>
class EmbeddedParallelTissueProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    EmbeddedParallelTissueProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                  std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Tissue")
    , couplingManager_(couplingManager)
    {}

    void enableCoupling(bool enable)
    { couplingEnabled_ = enable; }

    //! Dirichlet conditions on the whole boundary instead of no flow, which makes the coupled problem uniquely solvable
    void setDirichletBoundary(bool dirichlet)
    { dirichlet_ = dirichlet; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (dirichlet_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0 + 0.5*globalPos[0] - 0.25*globalPos[2]); }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = couplingManager().bulkPointSources(); }

    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume& scv) const
    {
        if constexpr (CouplingManager::couplingMode != Embedded1d3dCouplingMode::kernel)
        {
            const Scalar exchange = couplingEnabled_ ? embeddedParallelExchange(couplingManager(), source.id()) : 0.0;
            source = exchange*source.quadratureWeight()*source.integrationElement();
        }
        else
            source = 0.0;
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
        {
            if (!couplingEnabled_)
                return source;

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& ids = couplingManager().bulkSourceIds(eIdx, scv.indexInElement());
            const auto& weights = couplingManager().bulkSourceWeights(eIdx, scv.indexInElement());
            for (std::size_t i = 0; i < ids.size(); ++i)
                source[0] += embeddedParallelExchange(couplingManager(), ids[i])*weights[i];
            source[0] /= scv.volume()*elemVolVars[scv].extrusionFactor();
        }
        return source;
    }

    //! the flux scaling factor of the kernel coupling for a distance r, vessel radius R and kernel width rho
    Scalar fluxScalingFactor(const Scalar r, const Scalar R, const Scalar rho) const
    { return 1.0/(1.0 + 3.0*r + 5.0*R + 7.0*rho); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    bool couplingEnabled_ = true;
    bool dirichlet_ = false;
};

template<class TypeTag>
class EmbeddedParallelNetworkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    EmbeddedParallelNetworkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Network")
    , couplingManager_(couplingManager)
    {}

    void enableCoupling(bool enable)
    { couplingEnabled_ = enable; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = couplingManager().lowDimPointSources(); }

    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume& scv) const
    {
        const Scalar exchange = couplingEnabled_ ? embeddedParallelExchange(couplingManager(), source.id()) : 0.0;
        source = -exchange*source.quadratureWeight()*source.integrationElement();
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    bool couplingEnabled_ = true;
};

} // end namespace Dumux

#endif
