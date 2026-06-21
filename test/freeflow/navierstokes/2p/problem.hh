// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_TWOP_RISING_BUBBLE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_TWOP_RISING_BUBBLE_PROBLEM_HH

#include <algorithm>

#include <dumux/common/math.hh>
#include <cmath>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/vtkoutputmodule.hh>

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
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using DirichletValues = typename ParentType::DirichletValues;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::LocalView;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Vertex = typename GridGeometry::GridView::template Codim<dim>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class InitialShape { circle, square };

    // Surface-tension force discretization:
    //  - stress:       deviatoric Korteweg face flux γε(∇φ⊗∇φ)·n (well conditioned but
    //                  NOT balanced-force → steady spurious currents).
    //  - potential:    volumetric μ∇φ (balanced in the continuum but μ~γ/ε is large).
    //  - wellBalanced: volumetric -φ∇μ. Equivalent to μ∇φ up to a pure gradient that the
    //                  (modified) pressure absorbs. At CH equilibrium ∇μ≡0 element-wise, so
    //                  the discrete force is exactly zero → no spurious currents.
    enum class SurfaceTensionForm { stress, potential, wellBalanced };

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

        // Surface-tension force discretization. The balanced-force "wellBalanced" form
        // (-φ∇μ) gives a static bubble with (near) zero spurious currents; the "stress"
        // form is well conditioned but not balanced; "potential" is μ∇φ.
        // Back-compat: if SurfaceTensionForm is unset, fall back to UseKortewegStressForm.
        const auto stForm = getParam<std::string>("Problem.SurfaceTensionForm",
            getParam<bool>("Problem.UseKortewegStressForm", true) ? "stress" : "potential");
        if (stForm == "stress")
            surfaceTensionForm_ = SurfaceTensionForm::stress;
        else if (stForm == "potential")
            surfaceTensionForm_ = SurfaceTensionForm::potential;
        else if (stForm == "wellbalanced" || stForm == "wellBalanced")
            surfaceTensionForm_ = SurfaceTensionForm::wellBalanced;
        else
            DUNE_THROW(Dumux::ParameterException, "Unknown SurfaceTensionForm: " << stForm);

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
            // Hysing et al. (2009) FeatFlow rising-bubble benchmark boundary conditions:
            //   - no-slip (u = v = 0) on the top and bottom walls
            //   - free-slip on the left and right vertical walls
            //     (no penetration u_x = 0 + zero tangential shear stress)
            // This closes the box; with no through-flow the pressure is fixed only up to a
            // constant, which the optional bottom anchor below removes.
            const auto& pos = scv.dofPosition();
            if (isTop_(pos) || isBottom_(pos))
                values.setAllDirichlet(); // no-slip
            else
            {
                values.setDirichlet(Indices::velocityXIdx);      // no penetration
                values.setNeumann(Indices::momentumYBalanceIdx); // zero tangential shear
            }
        }
        else
        {
            values.setAllNeumann();
            // Optionally anchor the pressure level (the mass block has an all-zero pressure
            // column, so without a reference the saddle-point Jacobian can be singular).
            // The closed box leaves the pressure undetermined up to a constant.
            static const bool anchorPressure = getParam<bool>("Problem.AnchorPressureAtBottom", true);
            if (anchorPressure)
            {
                const auto& pos = scv.dofPosition();
                const bool isBottom = pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + 1e-6;
                if (isBottom)
                    values.setDirichlet(Indices::pressureIdx);
            }
        }

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
            // volumetric capillary force (used unless the stress face-flux form is active).
            if (surfaceTensionForm_ != SurfaceTensionForm::stress)
            {
                const auto ip = ipData(fvGeometry, scv);
                const auto grads = this->couplingManager().gradients(element, fvGeometry, ip);
                const auto values = this->couplingManager().values(element, fvGeometry, ip);
                if (surfaceTensionForm_ == SurfaceTensionForm::wellBalanced)
                {
                    // balanced-force form -φ∇μ: at CH equilibrium ∇μ≡0 → exactly zero force.
                    // The solver pressure is then the smooth modified pressure p* = p - μφ;
                    // the Laplace jump lives in μφ, not in the discrete pressure gradient.
                    const auto phi = values[CouplingManager::phaseFieldIdx];
                    source = -phi * grads[CouplingManager::chemicalPotentialIdx];
                }
                else // potential form μ∇φ
                {
                    const auto mu = values[CouplingManager::chemicalPotentialIdx];
                    source = mu * grads[CouplingManager::phaseFieldIdx];
                }
            }
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
        const Scalar phiIC = std::tanh((radius - dist) / (std::sqrt(2.0)*interfaceThickness_));
        initial[Indices::phaseFieldIdx] = phiIC;

        // Initialize the chemical potential to its φ-consistent equilibrium value, so the
        // coupled solve does NOT have to build μ from 0 in the first step (that transient
        // creates a huge ∇μ / μ and a spurious ~10^3 m/s velocity kick that displaces the
        // bubble and forces tiny time steps). For the tanh circle, inserting φ into
        // μ = (γ/ε)f'(φ) - γε∇²φ (with ∇² in polar coords) the double-well and normal
        // curvature terms cancel, leaving the geometric (Gibbs-Thomson) part:
        //     μ_eq(r) = γ (1 - φ²) / (√2 r),   γ = scaledSurfaceTension_, r = dist to centre.
        if (initialShape_ == InitialShape::circle && dist > 1e-8)
            initial[Indices::chemicalPotentialIdx] = scaledSurfaceTension_ * (1.0 - phiIC*phiIC) / (std::sqrt(2.0) * dist);
        else
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
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        const auto grads = this->couplingManager().gradients(element, fvGeometry, fluxVarCache.ipData());
        const auto& gradPhi = grads[CouplingManager::phaseFieldIdx];
        const auto normal = scvf.unitOuterNormal();
        using ResultType = std::decay_t<decltype(gradPhi)>;
        ResultType flux(0.0);

        // Korteweg capillary stress: T_cap n = γε (∇φ ⊗ ∇φ) n.
        // This is the well-conditioned surface-tension force (coefficient γε ~ O(σε)),
        // in contrast to the volumetric μ∇φ source (μ ~ γ/ε is large). Not balanced-force,
        // so it produces steady spurious currents; prefer SurfaceTensionForm = wellbalanced.
        if (surfaceTensionForm_ == SurfaceTensionForm::stress)
            flux += (scaledSurfaceTension_ * interfaceThickness_) * gradPhi * (gradPhi * normal);

        // Inertial conservative-form mass-flux correction:
        //   T_diff n = v (M grad(mu) . n) * (rho1 - rho2) / 2
        // Only meaningful with inertia enabled.
        if (enableInterfaceFlux_ && this->enableInertiaTerms())
        {
            const auto& shapeValues = fluxVarCache.shapeValues();
            const auto& gradChemicalPotential = grads[CouplingManager::chemicalPotentialIdx];
            ResultType v(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                v.axpy(shapeValues[localDof.index()][0], elemVolVars[localDof.index()].velocity());
            flux += this->mobility() * v * (gradChemicalPotential * normal) * this->mixtureDensityDerivative();
        }

        flux *= Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        return flux;
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
        Scalar bubbleVol = 0.0;     // ∫ (1+φ)/2 dV  (phase-1 / bubble volume)
        Scalar bubbleYmoment = 0.0; // ∫ (1+φ)/2 · y dV
        // Circularity (Hysing benchmark): c = perimeter of the area-equivalent circle / bubble
        // perimeter = 2*sqrt(pi*A)/P. Both A and P are taken from the SAME φ=0 isocontour (via
        // marching triangles) so the isoperimetric bound c<=1 holds (mixing the diffuse area
        // ∫(1+φ)/2 with the sharp φ=0 contour breaks it and gives spurious c>1). This is also
        // robust to the φ profile not fully saturating at ±1. On this LEFT-half domain the
        // contour is the left half of the interface, so the full-bubble perimeter and area are
        // each 2× the half-domain values.
        Scalar contourLenHalf = 0.0; // length of the φ=0 contour in the half domain
        Scalar areaHalfPos = 0.0;    // area where φ>0 (bubble) in the half domain, clipped at φ=0
        Scalar tvHalf = 0.0;         // ∫|∇φ| (total variation) - independent perimeter cross-check
        const auto& gg = this->gridGeometry();
        auto fvGeometry = localView(gg);
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bind(element);
            const auto numCorners = element.geometry().corners();
            std::vector<Scalar> phiC(numCorners);
            std::vector<GlobalPosition> posC(numCorners);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto phi = sol[scv.dofIndex()][Indices::phaseFieldIdx];
                const auto volume = scv.volume();
                totalMassLiquid += rho1_ * 0.5 * (1.0 + phi) * volume;
                totalMassGas    += rho2_ * 0.5 * (1.0 - phi) * volume;
                const auto bubbleFrac = 0.5 * (1.0 + phi) * volume;
                bubbleVol += bubbleFrac;
                bubbleYmoment += bubbleFrac * scv.dofPosition()[dimWorld-1];
                phiC[scv.localDofIndex()] = phi;
                posC[scv.localDofIndex()] = scv.dofPosition();
            }

            // Marching triangles: the φ=0 contour gives BOTH the interface length and the
            // enclosed (φ>0) area, clipped consistently at φ=0.
            if (numCorners == 3)
            {
                auto triArea = [](const GlobalPosition& a, const GlobalPosition& b, const GlobalPosition& c)
                { return 0.5*std::abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])); };
                auto edgeCross = [&](int a, int b)
                { const Scalar t = phiC[a]/(phiC[a]-phiC[b]); auto p = posC[a]; p.axpy(t, posC[b]-posC[a]); return p; };

                const int npos = int(phiC[0] > 0.0) + int(phiC[1] > 0.0) + int(phiC[2] > 0.0);
                const Scalar fullA = triArea(posC[0], posC[1], posC[2]);
                if (npos == 3)
                    areaHalfPos += fullA;
                else if (npos == 1) // one positive corner i: φ>0 region is the small corner triangle
                {
                    const int i = phiC[0] > 0.0 ? 0 : (phiC[1] > 0.0 ? 1 : 2);
                    const auto Xij = edgeCross(i, (i+1)%3);
                    const auto Xik = edgeCross(i, (i+2)%3);
                    areaHalfPos += triArea(posC[i], Xij, Xik);
                    contourLenHalf += (Xik - Xij).two_norm();
                }
                else if (npos == 2) // one negative corner k: φ>0 region is the complementary quad
                {
                    const int k = phiC[0] < 0.0 ? 0 : (phiC[1] < 0.0 ? 1 : 2);
                    const auto Xki = edgeCross(k, (k+1)%3);
                    const auto Xkj = edgeCross(k, (k+2)%3);
                    areaHalfPos += fullA - triArea(posC[k], Xki, Xkj);
                    contourLenHalf += (Xkj - Xki).two_norm();
                }

                // Independent perimeter cross-check: total variation ∫|∇φ| (also exact for the P1
                // field, but from the gradient magnitude, not the contour location). On the half
                // domain ∫|∇φ| equals the FULL-bubble perimeter IF φ is monotone across the
                // interface and saturates at ±1. A gap vs the contour perimeter therefore flags an
                // unsaturated/asymmetric diffuse interface, not a contour-extraction error.
                const Scalar twoA = (posC[1][0]-posC[0][0])*(posC[2][1]-posC[0][1])
                                  - (posC[2][0]-posC[0][0])*(posC[1][1]-posC[0][1]);
                if (std::abs(twoA) > 1e-30)
                {
                    const Scalar dphidx = ((phiC[1]-phiC[0])*(posC[2][1]-posC[0][1])
                                         - (phiC[2]-phiC[0])*(posC[1][1]-posC[0][1])) / twoA;
                    const Scalar dphidy = ((phiC[2]-phiC[0])*(posC[1][0]-posC[0][0])
                                         - (phiC[1]-phiC[0])*(posC[2][0]-posC[0][0])) / twoA;
                    tvHalf += std::hypot(dphidx, dphidy) * fullA;
                }
            }
        }
        const Scalar areaFull = 2.0 * areaHalfPos;
        const Scalar perimeterFull = 2.0 * contourLenHalf;
        const Scalar perimeterTV = tvHalf; // independent (total-variation) perimeter estimate
        const Scalar circularity = perimeterFull > 0.0 ? 2.0 * std::sqrt(M_PI * areaFull) / perimeterFull : 0.0;
        const Scalar circularityTV = perimeterTV > 0.0 ? 2.0 * std::sqrt(M_PI * areaFull) / perimeterTV : 0.0;

        std::cout << "\033[1;36m[mass] total = " << totalMassLiquid + totalMassGas
                  << " kg/m  |  liquid = " << totalMassLiquid
                  << " kg/m  |  gas = " << totalMassGas
                  << " kg/m\033[0m" << std::endl;
        std::cout << "\033[1;35m[bubble] centroid_y = " << bubbleYmoment / bubbleVol
                  << " m  |  volume = " << bubbleVol << " m^2/m\033[0m" << std::endl;
        std::cout << "\033[1;32m[bubble] circularity = " << circularity
                  << "  |  perimeter = " << perimeterFull << " m  |  area = " << areaFull << " m^2\033[0m" << std::endl;
        std::cout << "\033[1;32m[bubble] circularity_TV = " << circularityTV
                  << "  |  perimeter_TV = " << perimeterTV << " m  (independent cross-check)\033[0m" << std::endl;
    }

private:
    bool isTop_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - 1e-6; }
    bool isBottom_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + 1e-6; }
    Scalar energyScale_, mobility_, surfaceTension_, interfaceThickness_, scaledSurfaceTension_;
    Scalar rho1_, rho2_;
    Scalar eta1_, eta2_;
    Scalar g_;
    bool enableDampingForce_;
    Scalar dampingForceMagnitude_;
    bool enableInterfaceFlux_;
    SurfaceTensionForm surfaceTensionForm_;

    InitialShape initialShape_;
};

} // end namespace Dumux

#endif
