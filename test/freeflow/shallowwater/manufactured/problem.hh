// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief A manufactured solution of the shallow water equations for convergence studies
 */
#ifndef DUMUX_SHALLOWWATER_MANUFACTURED_PROBLEM_HH
#define DUMUX_SHALLOWWATER_MANUFACTURED_PROBLEM_HH

#include <string>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>

#include "exactsolution.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief A manufactured solution of the shallow water equations for convergence studies
 *
 * The flow of ShallowWaterTest::ManufacturedSolution on the unit square with the frequency `Problem.Frequency`
 * (default 0, steady), driven by its sources. At the boundary, the flux is the solution of the Riemann problem
 * between the inside state and the exact state at the face. The exact state is given the bed of the inside
 * state, with the water depth shifted such that the free surface elevation is the exact one. The initial
 * solution is the exact solution at the cell centers, scaled by one plus `Problem.InitialPerturbation` times a
 * smooth function of the position. Sources and boundary states are evaluated at the time set by setTime.
 */
template <class TypeTag>
class ManufacturedSolutionProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using NeumannFluxes = NumEqVector;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Solution = ShallowWaterTest::ManufacturedSolution<Scalar>;

public:
    ManufacturedSolutionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , name_(getParam<std::string>("Problem.Name"))
    , initialPerturbation_(getParam<Scalar>("Problem.InitialPerturbation", 0.0))
    , frequency_(getParam<Scalar>("Problem.Frequency", 0.0))
    {}

    const std::string& name() const
    { return name_; }

    //! Set the time at which sources and boundary states are evaluated
    void setTime(const Scalar t) const
    { time_ = t; }

    Scalar frequency() const
    { return frequency_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    template<class ElementFluxVariablesCache>
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        const auto inside = AdvectionType::insideFaceState(*this, fvGeometry, elemVolVars, scvf);
        const auto [h, u, v, z] = Solution::state(scvf.center(), time_, frequency_);
        const auto gravity = this->spatialParams().gravity(scvf.center());
        const auto riemannFlux = ShallowWater::riemannProblem(inside.waterDepth, h + z - inside.bedSurface,
                                                              inside.velocityX, u,
                                                              inside.velocityY, v,
                                                              inside.bedSurface, inside.bedSurface,
                                                              gravity, scvf.unitOuterNormal());
        NeumannFluxes values(0.0);
        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx] = riemannFlux[1];
        values[Indices::velocityYIdx] = riemannFlux[2];
        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        const auto s = Solution::source(globalPos, this->spatialParams().gravity(globalPos), time_, frequency_);
        return { s[0], s[1], s[2] };
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        using std::sin;
        const auto [h, u, v, z] = Solution::state(globalPos, time_, frequency_);
        const Scalar factor = 1.0 + initialPerturbation_*sin(3.0*globalPos[0] + 2.0*globalPos[1]);
        return { h*factor, u*factor, v*factor };
    }

private:
    std::string name_;
    Scalar initialPerturbation_;
    Scalar frequency_;
    mutable Scalar time_ = 0.0;
};

} // end namespace Dumux

#endif
