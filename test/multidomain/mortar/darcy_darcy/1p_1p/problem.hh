// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The problem for the single-phase darcy-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_DARCY_ONEP_TEST_PROBLEM_HH
#define DUMUX_MORTAR_DARCY_ONEP_TEST_PROBLEM_HH

#include <ranges>
#include <optional>
#include <unordered_map>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/facetgrid.hh>
#include <dumux/discretization/traceoperator.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/multidomain/mortar/subdomainproblembase.hh>

#include "property_declarations.hh"

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The problem for the single-phase darcy-darcy mortar-coupling test
 */
template<class TypeTag>
class DarcyProblem
: public PorousMediumFlowProblem<TypeTag>
, public Mortar::SubDomainFVProblemBase<
    GetPropType<TypeTag, Properties::GridGeometry>,
    GetPropType<TypeTag, Properties::MortarGrid>,
    GetPropType<TypeTag, Properties::MortarSolutionVector>
>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using MortarProblemBase = Mortar::SubDomainFVProblemBase<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::MortarGrid>,
        GetPropType<TypeTag, Properties::MortarSolutionVector>
    >;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>;

public:
    DarcyProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , MortarProblemBase(gridGeometry)
    {}

    template<typename GridVariables>
    auto assembleTraceVariables(std::size_t mortarId,
                                const GridVariables& gridVariables,
                                const SolutionVector& x) const
    {
        return this->getTraceOperator(mortarId).assembleScvfVariables(
            [&] (const auto& scvf, const auto& context) -> Dune::FieldVector<double, 1> {
                FluxVariables fluxVars;
                fluxVars.init(
                    *this,
                    context.fvGeometry().element(), context.fvGeometry(),
                    context.elemVolVars(), scvf, context.elemFluxVarsCache()
                );
                return {fluxVars.advectiveFlux(0, [] (const auto& vv) { return vv.mobility(0); })};
            },
            FVTraceOperatorContext<GridVariables, SolutionVector>{gridVariables, x}
        );
    }

    template<typename SubControlEntity>
    BoundaryTypes boundaryTypes(const Element& element, const SubControlEntity& sce) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        if (this->isOnMortarBoundary(element, sce))
            values.setAllDirichlet();
        return values;
    }

    template<typename SubControlEntity>
    PrimaryVariables dirichlet(const Element& element, const SubControlEntity& sce) const
    {
        if (this->isOnMortarBoundary(element, sce))
            return this->getTraceVariables(element, sce);
        if (useHomogeneousBoundaryConditions_)
            return PrimaryVariables(0);
        if constexpr (isCVFE)
            return PrimaryVariables(exactSolution_(sce.dofPosition()));
        else
            return PrimaryVariables(exactSolution_(sce.ipGlobal()));
    }

    void useHomogeneousBoundaryCondition(bool value)
    { useHomogeneousBoundaryConditions_ = value; }

private:
    double exactSolution_(const GlobalPosition& globalPos) const
    { return globalPos[0]*globalPos[1]; }

    bool useHomogeneousBoundaryConditions_ = false;
};

} // end namespace Dumux

#endif
