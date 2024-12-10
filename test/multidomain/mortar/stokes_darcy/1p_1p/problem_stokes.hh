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
 * \brief The problem for the stokes domain in the single-phase stokes-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_STOKES_DARCY_ONEP_TEST_STOKES_PROBLEM_HH
#define DUMUX_MORTAR_STOKES_DARCY_ONEP_TEST_STOKES_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/multidomain/mortar/subdomainproblembase.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The problem for the stokes domain in the single-phase stokes-darcy mortar-coupling test
 */
template<typename TypeTag, typename BaseProblem>
class StokesProblem
: public BaseProblem
, public Mortar::SubDomainFVProblemBase<
    GetPropType<TypeTag, Properties::GridGeometry>,
    GetPropType<TypeTag, Properties::MortarGrid>,
    GetPropType<TypeTag, Properties::MortarSolutionVector>
>
{
    using MortarProblemBase = Mortar::SubDomainFVProblemBase<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::MortarGrid>,
        GetPropType<TypeTag, Properties::MortarSolutionVector>
    >;

    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GG::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GG::LocalView;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using BoundaryTypes = typename BaseProblem::BoundaryTypes;
    using DirichletValues = typename BaseProblem::DirichletValues;
    using BoundaryFluxes = typename BaseProblem::BoundaryFluxes;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

public:
    using GridGeometry = GG;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    StokesProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                  std::shared_ptr<CouplingManager> couplingManager)
    : BaseProblem(gridGeometry, couplingManager)
    , MortarProblemBase{gridGeometry}
    {}

    template<typename GridVariables>
    auto assembleTraceVariables(std::size_t mortarId,
                                const GridVariables& gridVariables,
                                const SolutionVector& x) const
    {
        static_assert(BaseProblem::isMomentumProblem());
        return this->getTraceOperator(mortarId).assembleScvfVariables(
            [&] (const auto& scvf, const auto& context) -> Dune::FieldVector<double, 1> {
                return x[context.fvGeometry().scv(scvf.insideScvIdx()).dofIndex()]*scvf.directionSign();
            },
            FVTraceOperatorContext<GridVariables, SolutionVector>{gridVariables, x}
        );
    }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace scvf) const
    {
        BoundaryTypes values;
        if constexpr (BaseProblem::isMomentumProblem())
        {
            values.setAllDirichlet();
            if (onRightBoundary_(scvf.center()) or onLeftBoundary_(scvf.center()))
                values.setAllNeumann();
            else if (this->isOnMortarBoundary(element, scvf))
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();
        return values;
    }

    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    { return DirichletValues(0.0); }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (!BaseProblem::isMomentumProblem())
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf)
                                            *elemVolVars[scvf.insideScvIdx()].density()
                                            *scvf.unitOuterNormal();
        else
            if (this->isOnMortarBoundary(element, scvf))
                values[Indices::velocityYIdx] = this->getTraceVariables(element, scvf)[0]*scvf.directionSign();

        return values;
    }

    void useHomogeneousBoundaryCondition(bool value)
    { useHomogeneousBoundaryConditions_ = value; }

private:
    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + 1e-6; }

    bool onRightBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6; }

    bool onTopBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6; }

    bool onBottomBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6; }

    bool useHomogeneousBoundaryConditions_ = false;
};

} // end namespace Dumux

#endif
