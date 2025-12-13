// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_TWOP_RISING_BUBBLE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_TWOP_RISING_BUBBLE_PROBLEM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/grid/griddata.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/boundingboxtree.hh>

namespace Dumux {

/*!
 * \brief Test problem for the (Navier-) Stokes model
 *
 * Rising bubble in a resting fluid in a vertical column
 */
template <class TypeTag, class BaseProblem>
class ThreeDChannelTestProblem : public BaseProblem
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
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::LocalView;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Vertex = typename GridGeometry::GridView::template Codim<dim>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    ThreeDChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension", 24.5);
        interfaceThickness_ = getParam<Scalar>("Problem.InterfaceThickness", 0.08);

        // Mobility M controls diffusion rate in CH equation: flux = -M(φ²-1)²∇μ
        // Typical scaling: M ~ ε² for interface-controlled dynamics
        mobility_ = getParam<Scalar>("Problem.Mobility", 1.0);

        // Surface tension energy γ appears in: μ = -γε∇²φ + (γ/ε)f'(φ)
        // Standard CH uses γ = 3σ/(2√2) where σ is physical surface tension
        scaledSurfaceTension_ = surfaceTension_ * 3.0 / (2.0 * std::sqrt(2.0)) * getParam<Scalar>("Problem.SurfaceTensionFactor", 1.0);

        // Energy scale (γ/ε) for the double-well potential derivative: f'(φ) = (γ/ε)φ(φ²-1)
        energyScale_ = scaledSurfaceTension_ / interfaceThickness_ * getParam<Scalar>("Problem.EnergyScaleFactor", 1.0);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scv The sub control volume
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolume& scv) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (isTop_(scv.dofPosition()) || isBottom_(scv.dofPosition()))
                values.setAllDirichlet();
            else
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const

    {
        Sources source(0.0);
        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto& volVars = elemVolVars[scv];
            const auto phi = volVars.phaseField();
            // Double-well potential derivative: f'(φ) = (γ/ε)·φ(φ²-1)
            // f(φ) = (γ/ε)·(1/4)(φ²-1)² → f'(φ) = (γ/ε)·(1/2)·2φ(φ²-1)
            source[Indices::chemicalPotentialEqIdx] = volVars.chemicalPotential() - energyScale_*phi*(phi*phi - 1.0);
        }
        return source;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // no-flow/no-slip
        return DirichletValues(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
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
            if (isSide_(scvf.ipGlobal()))
            {
                values[Indices::conti0EqIdx] = 0.0;
                values[Indices::phaseFieldEqIdx] = 0.0;
            }
        }
        else
        {
            if (this->enableInertiaTerms())
            {
                if (isSide_(scvf.ipGlobal()))
                {
                    // advective term: vv*n
                    const auto& fluxVarCache = elemFluxVarsCache[scvf];
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto v = evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, scvf.ipGlobal());
                    const auto rho = this->couplingManager().density(element, fvGeometry, fluxVarCache.ipData());
                    values.axpy(rho*(v*scvf.unitOuterNormal()), v);

                    // stress tensor
                    Tensor gradV(0.0);
                    for (const auto& localDof : localDofs(fvGeometry))
                    {
                        const auto& volVars = elemVolVars[localDof];
                        for (int dir = 0; dir < dim; ++dir)
                            gradV[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(localDof.index()));
                    }

                    // get viscosity from the problem
                    const auto mu = this->couplingManager().effectiveViscosity(element, fvGeometry, fluxVarCache.ipData());
                    // compute -mu*gradV*n*dA
                    BoundaryFluxes momFlux = mv(gradV + getTransposed(gradV), scvf.unitOuterNormal());
                    momFlux *= -mu;

                    const auto pressure = this->couplingManager().pressure(element, fvGeometry, fluxVarCache.ipData());
                    momFlux += pressure * scvf.unitOuterNormal();

                    const auto normal = scvf.unitOuterNormal();
                    values += (momFlux * normal) * normal;
                }
            }
        }

        return values;
    }

    Scalar mobility() const
    {
        return mobility_;
    }

    Scalar surfaceTension() const
    {
        return scaledSurfaceTension_*interfaceThickness_;
    }

    /*!
     * \brief Evaluates the initial value for a control volume (momentum => velocity)
     */
    InitialValues initial(const SubControlVolume& scv) const
    {
        return InitialValues(0.0);
    }

    /*!
     * \brief Evaluates the initial value for a vertex (mass => pressure, phase field, chemical potential)
     */
    InitialValues initial(const Vertex& vertex) const
    {
        InitialValues initial(0.0);

        const auto pos = vertex.geometry().corner(0);

        const Scalar centerX = 0.5;
        const Scalar centerY = 0.5;

        const Scalar radius = 0.25;
        const Scalar dist = std::hypot(pos[0] - centerX, pos[1] - centerY);

        // Phase field: equilibrium tanh profile φ(r) = tanh((R-r)/ε)
        // For the double-well f(φ) = (γ/ε)·(1/4)(φ²-1)², the equilibrium requires:
        // -γε·∇²φ + (γ/ε)·φ(φ²-1) = μ
        initial[Indices::phaseFieldIdx] = std::tanh((radius - dist) / interfaceThickness_);
        initial[Indices::chemicalPotentialIdx] = 0.0;

        return initial;
    }

    GravityVector gravity() const
    { return {0.0, -0.98}; }

    /*!
     * \brief Returns the interface flux contribution due to the capillary forces
     */
    template<class Context>
    auto interfaceFlux(const Context& context) const
    {
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        const auto gradPhi = this->couplingManager().gradPhaseField(element, fvGeometry, fluxVarCache.ipData());
        const auto normal = scvf.unitOuterNormal();

        // T_cap = γε[∇φ⊗∇φ]
        return (scaledSurfaceTension_ * interfaceThickness_)
            * ((gradPhi * normal) * gradPhi)
            * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
    }

    Scalar mixtureViscosity(const Scalar phaseField) const
    {
        constexpr Scalar eta1 = 10.0;
        constexpr Scalar eta2 = 10.0;
        return 0.5 * ((1.0 + phaseField) * eta1 + (1.0 - phaseField) * eta2);
    }

    Scalar mixtureDensity(const Scalar phaseField) const
    {
        constexpr Scalar rho1 = 500.0;
        constexpr Scalar rho2 = 1000.0;
        return 0.5 * ((1.0 + phaseField) * rho1 + (1.0 - phaseField) * rho2);
    }

    Scalar mixtureDensityDerivative() const
    {
        constexpr Scalar rho1 = 500.0;
        constexpr Scalar rho2 = 1000.0;
        return 0.5 * (rho1  - rho2);
    }

private:
    bool isTop_(const GlobalPosition& pos, const Scalar eps = eps_) const
    { return pos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps; }

    bool isBottom_(const GlobalPosition& pos, const Scalar eps = eps_) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps; }

    bool isSide_(const GlobalPosition& pos) const
    { return pos[0] < this->gridGeometry().bBoxMin()[0] + eps_ ||
             pos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_ = 1e-6;
    Scalar energyScale_, mobility_, surfaceTension_, interfaceThickness_, scaledSurfaceTension_;
};

} // end namespace Dumux

#endif
