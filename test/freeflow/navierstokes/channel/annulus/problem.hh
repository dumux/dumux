// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Radial flow from a line source between two coaxial cylinders
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ANNULUS_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ANNULUS_PROBLEM_HH

#include <bitset>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Radial flow from a line source between two coaxial cylinders
 *
 * The velocity \f$ u_r = A/r \f$ is divergence-free and irrotational, so its viscous force vanishes
 * and the pressure follows from the inertia alone, \f$ p = p_0 - \rho A^2/(2 r^2) \f$. Both cylinders
 * prescribe the velocity, the planes normal to the axis are planes of symmetry.
 */
template <class TypeTag, class BaseProblem>
class AnnulusTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using InitialValues = typename ParentType::InitialValues;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename ModelTraits::Indices;

    AnnulusTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        sourceStrength_ = getParam<Scalar>("Problem.SourceStrength");
        referencePressure_ = getParam<Scalar>("Problem.ReferencePressure");
        density_ = getParam<Scalar>("Component.LiquidDensity");
        enableInertiaTerms_ = getParam<bool>("Problem.EnableInertiaTerms");
        eps_ = 1e-7*(this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0]);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (onCylinder_(globalPos))
                values.setAllDirichlet();
            else
            {
                values.setDirichlet(Indices::velocityYIdx);
                values.setNeumann(Indices::momentumXBalanceIdx);
            }
        }
        else
            values.setAllNeumann();

        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
            if (onCylinder_(scvf.ipGlobal()))
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }

    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        InitialValues values(0.0);
        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::pressureIdx] = referencePressure_;
        return values;
    }

    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values(0.0);
        const auto r = globalPos[0];

        if constexpr (ParentType::isMomentumProblem())
            values[Indices::velocityXIdx] = sourceStrength_/r;
        else
        {
            values[Indices::pressureIdx] = referencePressure_;
            if (enableInertiaTerms_)
                values[Indices::pressureIdx] -= 0.5*density_*sourceStrength_*sourceStrength_/(r*r);
        }

        return values;
    }

    //! The pressure is only defined up to a constant, it is fixed in one cell
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    std::bitset<DirichletValues::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<DirichletValues::dimension> values;
        if (scv.dofIndex() == 0)
            values.set(Indices::pressureIdx);
        return values;
    }

    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return analyticalSolution(scv.dofPosition()); }

private:
    bool onCylinder_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_
            || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    Scalar sourceStrength_;
    Scalar referencePressure_;
    Scalar density_;
    bool enableInertiaTerms_;
    Scalar eps_;
};

} // end namespace Dumux

#endif
