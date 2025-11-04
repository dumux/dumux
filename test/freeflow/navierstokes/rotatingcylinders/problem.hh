// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH

#include <cmath>
#include <numeric>

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/grid/griddata.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \brief Test problem for the (Navier-) Stokes model in a 3D channel
 *
 * Benchmark case from
 *   Turek, Schaefer et a (1996) Benchmark computations of laminar flow around cylinder.
 *   Flow Simulation with High-Performance Computers II,
 *   Notes on Numerical Fluid Mechanics 52, 547-566, Vieweg 1996
 *   https://doi.org/10.1007/978-3-322-89849-4_39
 */
template <class TypeTag, class BaseProblem>
class DFGChannelTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using DirichletValues = typename ParentType::DirichletValues;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    DFGChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        density_ = getParam<Scalar>("Component.LiquidDensity");
        viscosity_ = getParam<Scalar>("Component.LiquidDynamicViscosity");
        const auto radii = getParam<std::vector<Scalar>>("Grid.Radial0");
        radius1_ = radii.front();
        radius2_ = radii.back();
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet();
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary face
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary face
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf)
                                           *density_*scvf.unitOuterNormal();

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     * \param time The current simulation time
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos) const
    {
        const Scalar omega1 = 1e2; // [rad/s]
        const Scalar omega2 = 0;
        const auto mu = omega2/omega1;
        const auto r2byr1 = radius2_*radius2_/(radius1_*radius1_);
        const auto a = omega1*(1 - mu*r2byr1)/(1 - r2byr1);
        const auto b = radius1_*radius1_*omega1*(1 - mu)/(1 - 1/r2byr1);

        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::sqrt;
        const auto radius = sqrt(x*x + y*y);

        DirichletValues values;
        if constexpr (ParentType::isMomentumProblem())
        {
            const auto uRad = a*radius + b/radius;
            values[Indices::velocityXIdx] = -y*uRad/radius;
            values[Indices::velocityYIdx] = x*uRad/radius;
        }
        else
        {
            using std::log;
            values[Indices::pressureIdx] = 0.5*a*a*radius*radius + 2*a*b*log(radius) + 0.5*b*b/(radius*radius);
        }

        return values;
    }

private:
    Scalar radius1_, radius2_, density_, viscosity_;
};

} // end namespace Dumux

#endif
