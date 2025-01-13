// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

template<class GridView>
class MeshMotion
{
    using GridGeometry = BoxFVGridGeometry<double, GridView, true>;
    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
public:

    MeshMotion(const GridGeometry& gridGeometry)
    : displacement_(gridGeometry.numDofs()),
      displacementGridGeometry_(gridGeometry)
    {
        displacement_ = 0.0;
    }

    void deform()
    {
        for (const auto& vertex : vertices(displacementGridGeometry_.gridView()))
        {
            const auto vIdx = displacementGridGeometry_.dofMapper().index(vertex);
            const auto globalPos = vertex.geometry().corner(0);
            displacement_[vIdx] = displacement(globalPos);
        }
    }

    const auto displacement(const GlobalPosition& globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        const auto dx = 0.0;
        const auto dy = [&]{
            if (y < 0.1)
                return -0.002*(1.0/(y - 0.2) + 5.0);
            else if (y > 0.1 && y <= 0.3)
                return -0.01*(y - 0.2)/0.1;
            else
                return -0.002*(1.0/(y - 0.2) - 4.7619047619);
        }();

        const auto smooth = x < 0.5 ? 0.0 : (x > 1.5 ? 0.0 : std::sin(M_PI*(x - 0.5)));

        return GlobalPosition{dx, 5*dy*smooth};
    }

    template<class FVGeometry>
    double detF(const FVGeometry& fvGeometry, const FVGeometry::SubControlVolume& scv) const
    {
        return detF(fvGeometry, scv.center());
    }

    template<class FVGeometry>
    double detF(const FVGeometry& fvGeometry, const FVGeometry::SubControlVolumeFace& scvf) const
    {
        return detF(fvGeometry, scvf.ipGlobal());
    }

    template<class FVGeometry>
    double detF(const FVGeometry& fvGeometry, const GlobalPosition& globalPos) const
    {
        return F(fvGeometry, globalPos).determinant();
    }

    template<class FVGeometry>
    Dune::FieldMatrix<double, dimWorld, dimWorld> F(const FVGeometry& fvGeometry, const FVGeometry::SubControlVolumeFace& scvf) const
    {
        return F(fvGeometry, scvf.ipGlobal());
    }

    template<class FVGeometry>
    Dune::FieldMatrix<double, dimWorld, dimWorld> F(const FVGeometry& fvGeometry, const GlobalPosition& globalPos) const
    {
        const auto elemSol = elementSolution(fvGeometry.element(), displacement_, displacementGridGeometry_);
        const auto& element = fvGeometry.element();

        const auto F = [&,this]{
            Dune::FieldMatrix<double, dimWorld, dimWorld> F(0.0);
            const auto gradU = evalGradients(element, element.geometry(), displacementGridGeometry_, elemSol, globalPos);
            for (int i = 0; i < dimWorld; ++i)
            {
                F[i][i] += 1.0;
                for (int j = 0; j < dimWorld; ++j)
                    F[i][j] += gradU[i][j];
            }
            return F;
        }();

        return F;
    }

    const Dune::BlockVector<GlobalPosition>& displacement() const
    { return displacement_; }

private:
    Dune::BlockVector<GlobalPosition> displacement_;
    const GridGeometry& displacementGridGeometry_;
};

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
    using GridView = typename GridGeometry::GridView;
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
    DFGChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager, std::shared_ptr<MeshMotion<GridView>> meshMotion)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    , meshMotion_(meshMotion)
    {
        density_ = getParam<Scalar>("Component.LiquidDensity");
        viscosity_ = getParam<Scalar>("Component.LiquidDynamicViscosity");

        for (int i = 0; i < dimWorld; ++i)
            domainSize_[i] = gridGeometry->bBoxMax()[i] -  gridGeometry->bBoxMin()[i];
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
        {
            if (isOutlet_(globalPos))
                values.setAllNeumann();
            else
                values.setAllDirichlet();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary face
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            return inflowVelocity_(globalPos[1]);
        else
            return DirichletValues(0.0);
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
        {
            values[Indices::conti0EqIdx]
                = this->faceVelocity(element, fvGeometry, scvf)
                    * density_ * scvf.unitOuterNormal();
        }

        return values;
    }

    template<class ScvOrScvf>
    Scalar density(const Element&,
                   const FVElementGeometry&,
                   const ScvOrScvf&,
                   const bool isPreviousTimeStep = false) const
    {
        return density_;
    }

    template<class ScvOrScvf>
    Scalar effectiveViscosity(const Element&,
                              const FVElementGeometry&,
                              const ScvOrScvf&) const
    {
        return viscosity_;
    }

    const MeshMotion<GridView>& meshMotion() const
    { return *meshMotion_; }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onCylinderBoundary_(const GlobalPosition& globalPos) const
    { return std::hypot(globalPos[0] - 0.2, globalPos[1] - 0.2) < 0.06; }

    DirichletValues inflowVelocity_(Scalar y) const
    {
        constexpr Scalar maxVelocity = 0.3;
        return { 4*maxVelocity*y*(domainSize_[1]-y)/(domainSize_[1]*domainSize_[1]), 0.0 };
    }

    static constexpr Scalar eps_ = 1e-10;

    std::shared_ptr<MeshMotion<GridView>> meshMotion_;

    Scalar density_, viscosity_;
    std::array<Scalar, dimWorld> domainSize_;
};

} // end namespace Dumux

#endif
