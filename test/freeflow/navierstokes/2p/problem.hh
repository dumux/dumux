// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_TWOP_RISING_BUBBLE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_TWOP_RISING_BUBBLE_PROBLEM_HH

#include <dune/common/float_cmp.hh>

#include <algorithm>

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

    enum class InitialShape { circle, square };

public:
    ThreeDChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension", 24.5);
        interfaceThickness_ = getParam<Scalar>("Problem.InterfaceThickness", 0.08);

        // Mobility M controls diffusion rate in CH equation: flux = -M∇μ
        // Typical scaling: M ~ ε² for interface-controlled dynamics
        mobility_ = getParam<Scalar>("Problem.Mobility", 1.0);

        // Surface tension energy γ appears in: μ = -γε∇²φ + (γ/ε)f'(φ)
        // For f(φ) = (1/4)(φ²-1)², the surface energy integral gives σ_eff = γ·(2√2/3),
        // so to recover the physical surface tension σ, we need γ = 3σ/(2√2).
        scaledSurfaceTension_ = surfaceTension_ * 3.0 / (2.0 * std::sqrt(2.0)) * getParam<Scalar>("Problem.SurfaceTensionFactor", 1.0);

        // Energy scale (γ/ε) for the double-well potential derivative: f'(φ) = (γ/ε)φ(φ²-1)
        energyScale_ = scaledSurfaceTension_ / interfaceThickness_ * getParam<Scalar>("Problem.EnergyScaleFactor", 1.0);

        // fluid parameters
        rho1_ = getParam<Scalar>("Problem.Density1", 1.0);
        rho2_ = getParam<Scalar>("Problem.Density2", 1000.0);
        eta1_ = getParam<Scalar>("Problem.Viscosity1", 1.0);
        eta2_ = getParam<Scalar>("Problem.Viscosity2", 10.0);
        g_ = getParam<Scalar>("Problem.Gravity", 0.98);

        enableDampingForce_ = getParam<bool>("Problem.EnableDampingForce", false);
        dampingForceMagnitude_ = getParam<Scalar>("Problem.DampingForceMagnitude", 1e5);
        enableInterfaceFlux_ = getParam<bool>("Problem.EnableInterfaceFlux", true);

        const auto shapeStr = getParam<std::string>("Problem.InitialShape", "circle");
        if (shapeStr == "circle")
            initialShape_ = InitialShape::circle;
        else if (shapeStr == "square")
            initialShape_ = InitialShape::square;
        else
            DUNE_THROW(Dumux::ParameterException, "Unknown InitialShape: " << shapeStr);
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
            if (isTop_(scv.dofPosition()))
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
        else
        {
            const auto ip = ipData(fvGeometry, scv);
            const auto grads = this->couplingManager().gradients(element, fvGeometry, ip);
            const auto mu = this->couplingManager().values(element, fvGeometry, ip)[CouplingManager::chemicalPotentialIdx];
            // volumetric Korteweg force (chemical-potential / phase-field gradient form)
            source = mu * grads[CouplingManager::phaseFieldIdx];
        }

        if constexpr (ParentType::isMomentumProblem())
        {
            // optional sponge-layer damping (off by default)
            if (enableDampingForce_ && scv.dofPosition()[1] < 0.0)
            {
                const auto v = elemVolVars[scv].velocity();
                source -= dampingForceMagnitude_ * v.two_norm() * v;
            }

            // Boussinesq buoyancy reference: combined with the base class' +rho_mix * g
            // this yields the effective body force (rho_mix - rho2) * g.
            source -= rho2_ * this->gravity();
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

            else if (isBottom_(scvf.ipGlobal()))
            {
                values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf)*scvf.unitOuterNormal();
                values[Indices::phaseFieldEqIdx] = 0.0;
            }
        }
        else
        {
            if (isSide_(scvf.ipGlobal()))
            {
                const auto& fluxVarCache = elemFluxVarsCache[scvf];
                if (this->enableInertiaTerms())
                {
                    // advective term: vv*n
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto v = evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, scvf.ipGlobal());
                    const auto rho = this->couplingManager().density(element, fvGeometry, fluxVarCache.ipData());
                    values.axpy(rho*(v*scvf.unitOuterNormal()), v);
                }

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

            else if (isBottom_(scvf.ipGlobal()))
            {
                const auto& fluxVarCache = elemFluxVarsCache[scvf];
                if (this->enableInertiaTerms())
                {
                    // advective term: vv*n
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto v = evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, scvf.ipGlobal());
                    const auto rho = this->couplingManager().density(element, fvGeometry, fluxVarCache.ipData());
                    values.axpy(rho*(v*scvf.unitOuterNormal()), v);
                }

                // zero pressure at bottom
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
        // Returns γε, the coefficient in front of the Laplacian in the chem. potential
        return scaledSurfaceTension_ * interfaceThickness_;
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

        const Scalar dist = initialShape_ == InitialShape::circle ? std::hypot(pos[0] - centerX, pos[1] - centerY)
            : std::max(std::abs((pos[0] - centerX)), std::abs((pos[1] - centerY)));

        // The tanh provides a smooth
        // diffuse-interface transition of thickness ~interfaceThickness_.
        initial[Indices::phaseFieldIdx] = std::tanh((radius - dist) / (std::sqrt(2.0)*interfaceThickness_));
        initial[Indices::chemicalPotentialIdx] = 0.0;

        return initial;
    }

    GravityVector gravity() const
    { return {0.0, -g_}; }

    /*!
     * \brief Returns the interface flux contribution due to the capillary forces
     */
    template<class Context>
    auto interfaceFlux(const Context& context) const
    {
        const auto& fvGeometry = context.fvGeometry();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];
        using ResultType = std::decay_t<decltype(this->couplingManager().gradients(
            context.element(), fvGeometry, fluxVarCache.ipData())[CouplingManager::phaseFieldIdx])>;

        // Korteweg / mass-flux interface contribution.
        // Off by default for pure Stokes runs; the volumetric Korteweg source
        // (mu * grad(phi)) already represents the capillary force.
        // Only the inertial conservative-form correction
        //   T_diff n = v (M grad(mu) . n) * (rho1 - rho2) / 2
        // is implemented here. Requires inertia to be physically meaningful.
        if (!enableInterfaceFlux_ || !this->enableInertiaTerms())
            return ResultType(0.0);

        const auto& element = context.element();
        const auto& elemVolVars = context.elemVolVars();
        const auto& shapeValues = fluxVarCache.shapeValues();

        const auto grads = this->couplingManager().gradients(element, fvGeometry, fluxVarCache.ipData());
        const auto& gradChemicalPotential = grads[CouplingManager::chemicalPotentialIdx];
        const auto normal = scvf.unitOuterNormal();
        ResultType v(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            v.axpy(shapeValues[localDof.index()][0], elemVolVars[localDof.index()].velocity());

        return this->mobility() * v * (gradChemicalPotential * normal) * this->mixtureDensityDerivative()
             * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
    }

    Scalar mixtureViscosity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, -1.0, 1.0);
        return 0.5 * ((1.0 + phi) * eta1_ + (1.0 - phi) * eta2_);
    }

    Scalar mixtureDensity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, -1.0, 1.0);
        return 0.5 * ((1.0 + phi) * rho1_ + (1.0 - phi) * rho2_);
    }

    Scalar mixtureDensityDerivative() const
    {
        return 0.5 * (rho1_  - rho2_);
    }

    template<class Solution, class GridVariables>
    void printMassBalanceSummary(const Solution& sol, const GridVariables& gridVariables) const
    {
        Scalar totalMassLiquid = 0.0;
        Scalar totalMassGas = 0.0;
        const auto& gg = this->gridGeometry();
        auto fvGeometry = localView(gg);
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto phi = sol[scv.dofIndex()][Indices::phaseFieldIdx];
                const auto volume = scv.volume();
                totalMassLiquid += rho1_ * 0.5 * (1.0 + phi) * volume;
                totalMassGas    += rho2_ * 0.5 * (1.0 - phi) * volume;
            }
        }

        std::cout << "\033[1;36m[mass] total = " << totalMassLiquid + totalMassGas
                  << " kg/m  |  liquid = " << totalMassLiquid
                  << " kg/m  |  gas = " << totalMassGas
                  << " kg/m\033[0m" << std::endl;
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
    Scalar rho1_, rho2_;
    Scalar eta1_, eta2_;
    Scalar g_;
    bool enableDampingForce_;
    Scalar dampingForceMagnitude_;
    bool enableInterfaceFlux_;

    InitialShape initialShape_;
};

} // end namespace Dumux

#endif
