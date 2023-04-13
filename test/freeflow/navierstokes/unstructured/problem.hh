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
        if constexpr (ParentType::isMomentumProblem())
        {
            if (this->enableInertiaTerms())
            {
                if (isOutlet_(scvf.ipGlobal()))
                {
                    // advective term: vv*n
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto v = evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, scvf.ipGlobal());
                    values.axpy(density_*(v*scvf.unitOuterNormal()), v);
                }
            }
        }
        else
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

    //! Computes pressure difference benchmark indicator
    //! (only works for correct domain size 2.2 x 0.41 and boundary conditions, Re=20)
    template<class SolutionVector, class GridVariables, class P = ParentType, typename std::enable_if_t<!P::isMomentumProblem(), int> = 0>
    Scalar evalPressureDifference(const GridVariables& gridVariables, const SolutionVector& p) const
    {
        const auto& gg = gridVariables.gridGeometry();
        const auto& tree = gg.boundingBoxTree();
        GlobalPosition evalPoint1({0.15, 0.2});
        GlobalPosition evalPoint2({0.25, 0.2});
        const Scalar stepSize = 0.001;

        const auto evalPressure = [&](auto pos, bool backwards)
        {
            auto entities = intersectingEntities(pos, tree);
            while (entities.empty())
            {
                pos[0] += backwards ? -stepSize : stepSize;
                entities = intersectingEntities(pos, tree);
            }

            const auto element = gg.element(entities[0]);
            const auto fvGeometry = localView(gg).bindElement(element);
            const auto elemVolVars = localView(gridVariables.curGridVolVars()).bindElement(element, fvGeometry, p);
            const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
            return evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, pos);
        };

        return evalPressure(evalPoint1, true) - evalPressure(evalPoint2, false);
    }

    //! Computes drag and lift coefficient benchmark indicator
    //! (only works for correct domain size 2.2 x 0.41 and boundary conditions, Re=20)
    template<class SolutionVector, class GridVariables, class P = ParentType, typename std::enable_if_t<P::isMomentumProblem(), int> = 0>
    auto evalDragAndLiftCoefficient(const GridVariables& gridVariables, const SolutionVector& v) const
    {
        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVariables.curGridVolVars());

        BoundaryFluxes forces(0.0);
        Scalar cylinderSurface = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            if (fvGeometry.hasBoundaryScvf())
            {
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (scvf.boundary() && onCylinderBoundary_(scvf.ipGlobal()))
                    {
                        elemVolVars.bind(element, fvGeometry, v);
                        Dune::FieldVector<GlobalPosition, 2> stress(GlobalPosition(0.0));

                        const auto geometry = element.geometry();
                        const auto isgeometry = fvGeometry.geometry(scvf);
                        const auto& quad = Dune::QuadratureRules<Scalar, std::decay_t<decltype(isgeometry)>::mydimension>::rule(isgeometry.type(), 5);
                        for (auto&& qp : quad)
                        {
                            const auto ipGlobal = isgeometry.global(qp.position());
                            const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                            auto stressLocal = evalGradients(element, geometry, fvGeometry.gridGeometry(), elemSol, ipGlobal);
                            stressLocal *= qp.weight()*isgeometry.integrationElement(qp.position());
                            stress = stress + stressLocal;
                        }
                        stress *= viscosity_;
                        stress /= scvf.area();

                        BoundaryFluxes normalStress(0.0);
                        for (int dir = 0; dir < dim; ++dir)
                        {
                            // Add pressure contribuition
                            stress[dir][dir] -= this->pressure(element, fvGeometry, scvf);
                            normalStress[dir] = stress[dir]*scvf.unitOuterNormal();
                        }
                        forces.axpy(scvf.area(), normalStress);

                        cylinderSurface += scvf.area();
                    }
                }
            }
        }

        const Scalar dragForce = -forces[0];
        const Scalar liftForce = -forces[1];

        const Scalar uMean = 0.2;
        const Scalar lChar = 0.1;
        const Scalar coefficientFactor = 2.0/(uMean*uMean*lChar);

        // sanity check
        std::cout << "Reynolds number: " << uMean*lChar/viscosity_ << std::endl;
        std::cout << "Cylinder surface: " << cylinderSurface
                  << " (reference: 0.31415926535)"
                  << std::endl;

        return std::make_pair( coefficientFactor*dragForce, coefficientFactor*liftForce );
    }

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
    Scalar density_, viscosity_;
    std::array<Scalar, dimWorld> domainSize_;
};

} // end namespace Dumux

#endif
