// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief A simple 1p problem on a 1d network embedded in a bulk domain, coupled by point sources.
 *        The spatial parameters (carrying the per-element radius) are injected, so the radius can
 *        come from a DGF file and be migrated onto the distributed grid.
 */
#ifndef DUMUX_NETWORK_PROBLEM_HH
#define DUMUX_NETWORK_PROBLEM_HH

#include <cmath>
#include <memory>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief 1p flow on a 1d network, coupled to the bulk by point sources.
 */
template <class TypeTag>
class NetworkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    template<class SpatialParams>
    NetworkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<SpatialParams> spatialParams,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, spatialParams, "Vessel")
    , couplingManager_(couplingManager)
    {
        name_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    const std::string& name() const
    { return name_; }

    //! Dirichlet at the network's outer (boundary) dofs; a depth-varying field drives flow.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 1.0 + globalPos[2];
        return values;
    }

    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    //! solution-dependent coupling source (Peaceman-type exchange term)
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume& scv) const
    {
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];
        const Scalar radius = this->couplingManager().radius(source.id());
        const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
        const Scalar sourceValue = beta*(pressure3D - pressure1D);
        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::string name_;
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
