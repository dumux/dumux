// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Rising-bubble Cahn-Hilliard / Navier-Stokes problem for the P2-P1-P2-P2 Taylor-Hood scheme
 *        (experimental "new" CVFE interface: constraint-based Dirichlet BCs).
 */
#ifndef DUMUX_TEST_FF_NAVIERSTOKES_2P_P2_PROBLEM_HH
#define DUMUX_TEST_FF_NAVIERSTOKES_2P_P2_PROBLEM_HH

#include <algorithm>
#include <cmath>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

template <class TypeTag, class BaseProblem>
class RisingBubbleP2Problem : public BaseProblem
{
    using ParentType = BaseProblem;

public:
    // Surface-tension force discretization (matches the pq1bubble baseline's
    // ThreeDChannelTestProblem::SurfaceTensionForm one folder up):
    //  - stress:       Korteweg face/quadrature flux gamma*eps*(gradC x gradC).n (well
    //                   conditioned but NOT balanced-force -> steady spurious currents).
    //  - potential:    volumetric mu*gradC (balanced in the continuum but mu~gamma/eps is large).
    //  - wellBalanced: volumetric -c*gradMu. Equivalent to mu*gradC up to a pure gradient
    //                  absorbed by the pressure; at CH equilibrium gradMu=0 elementwise, so the
    //                  discrete force is exactly zero -> no spurious currents. This is what the
    //                  validated Hysing benchmark run (../params.input) uses.
    enum class SurfaceTensionForm { stress, potential, wellBalanced };

private:

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using DirichletValues = typename ParentType::DirichletValues;
    using InitialValues = typename ParentType::InitialValues;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<typename GridGeometry::GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    RisingBubbleP2Problem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension", 24.5);
        interfaceThickness_ = getParam<Scalar>("Problem.InterfaceThickness", 0.005);
        mobility_ = getParam<Scalar>("Problem.Mobility", 5e-6);
        scaledSurfaceTension_ = surfaceTension_ * 3.0 / (2.0 * std::sqrt(2.0));
        energyScale_ = scaledSurfaceTension_ / interfaceThickness_;
        rho1_ = getParam<Scalar>("Problem.Density1", 100.0);
        rho2_ = getParam<Scalar>("Problem.Density2", 1000.0);
        eta1_ = getParam<Scalar>("Problem.Viscosity1", 1.0);
        eta2_ = getParam<Scalar>("Problem.Viscosity2", 10.0);
        g_ = getParam<Scalar>("Problem.Gravity", 0.98);
        radius_ = getParam<Scalar>("Problem.BubbleRadius", 0.25);

        const auto stForm = getParam<std::string>("Problem.SurfaceTensionForm", "wellbalanced");
        if (stForm == "stress")
            surfaceTensionForm_ = SurfaceTensionForm::stress;
        else if (stForm == "potential")
            surfaceTensionForm_ = SurfaceTensionForm::potential;
        else if (stForm == "wellbalanced" || stForm == "wellBalanced")
            surfaceTensionForm_ = SurfaceTensionForm::wellBalanced;
        else
            DUNE_THROW(Dumux::ParameterException, "Unknown SurfaceTensionForm: " << stForm);

        appendConstraints_();
    }

    //! all boundary equations are flux (natural / zero-flux); Dirichlet is enforced via constraints()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if constexpr (ParentType::isMomentumProblem())
            values.setAllFluxBoundary();
        else
            values.setFluxBoundary(Indices::conti0EqIdx);
        return values;
    }

    //! closed box: zero boundary flux everywhere (no through-flow)
    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    { return BoundaryFluxes(0.0); }

    template<class ElementVariables, class ElementFluxVariablesCache, class FaceIpData>
    BoundaryFluxes boundaryFlux(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const ElementFluxVariablesCache& elemFluxVarsCache,
                                const FaceIpData& faceIpData) const
    { return boundaryFlux(fvGeometry, elemVars, faceIpData); }

    //! Dirichlet velocity (used as constraint values); zero (no-slip / no-penetration)
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return DirichletValues(0.0); }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    //! initial condition (position based -> both vertex and P2 edge dofs)
    InitialValues initialAtPos(const GlobalPosition& pos) const
    {
        InitialValues values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar cx = 0.5, cy = 0.5;
            const Scalar dist = std::hypot(pos[0] - cx, pos[1] - cy);
            const Scalar phiIC = std::tanh((radius_ - dist)/(std::sqrt(2.0)*interfaceThickness_));
            values[Indices::phaseFieldIdx] = phiIC;
            values[Indices::chemicalPotentialIdx] = analyticChemicalPotential(pos);
        }
        return values; // mass: pressure = 0
    }

    //! The analytic Gibbs-Thomson equilibrium chemical potential of the initial tanh circle,
    //! mu_eq(r) = gamma~ (1 - phi^2) / (sqrt(2) r). Used for the IC and (diagnostically) to
    //! prescribe mu directly instead of solving the mu-equation (see prescribeChemicalPotential).
    Scalar analyticChemicalPotential(const GlobalPosition& pos) const
    {
        const Scalar cx = 0.5, cy = 0.5;
        const Scalar dist = std::hypot(pos[0] - cx, pos[1] - cy);
        const Scalar phi = std::tanh((radius_ - dist)/(std::sqrt(2.0)*interfaceThickness_));
        return (dist > 1e-8) ? scaledSurfaceTension_*(1.0 - phi*phi)/(std::sqrt(2.0)*dist) : 0.0;
    }

    //! Diagnostic switch: when true, the mu-definition equation is replaced by a pin
    //! mu = analyticChemicalPotential(x), so the momentum capillary force and the CH mobility
    //! flux both see a smooth, correct mu. Separates "mu-solve is producing a bad mu" from
    //! "the velocity/pressure response to a correct mu is wrong". Default off.
    bool prescribeChemicalPotential() const
    {
        static const bool p = getParam<bool>("Problem.PrescribeChemicalPotential", false);
        return p;
    }

    // ---- material / CH accessors used by the momentum-CHNS local residual ----
    Scalar mobility() const { return mobility_; }
    Scalar surfaceTension() const { return scaledSurfaceTension_ * interfaceThickness_; } // sigma~ eps
    Scalar energyScale() const { return energyScale_; }                                   // sigma~/eps
    Scalar referenceDensity() const { return rho2_; }
    SurfaceTensionForm surfaceTensionForm() const { return surfaceTensionForm_; }
    Scalar mixtureDensity(const Scalar phaseField) const
    { const auto p = std::clamp(phaseField, -1.0, 1.0); return 0.5*((1.0+p)*rho1_ + (1.0-p)*rho2_); }
    //! d rho / d c of the (unclamped) linear blend, for the Abels-Garcke-Grün momentum flux
    Scalar mixtureDensityDerivative() const { return 0.5*(rho1_ - rho2_); }
    Scalar mixtureViscosity(const Scalar phaseField) const
    { const auto p = std::clamp(phaseField, -1.0, 1.0); return 0.5*((1.0+p)*eta1_ + (1.0-p)*eta2_); }
    GravityVector gravity() const { return {0.0, -g_}; }

private:
    bool isTop_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - 1e-6; }
    bool isBottom_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + 1e-6; }

    //! Build the Dirichlet constraints: no-slip velocity on top/bottom, no-penetration velocityX on
    //! the free-slip side walls, and (mass) a single pressure anchor at the bottom-left corner.
    void appendConstraints_()
    {
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            if constexpr (ParentType::isMomentumProblem())
            {
                for (const auto& boundaryFace : boundaryFaces(fvGeometry))
                {
                    for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                    {
                        const auto& pos = ipData(fvGeometry, localDof).global();
                        ConstraintInfo info;
                        info.set(Indices::velocityXIdx);         // no penetration (all walls)
                        if (isTop_(pos) || isBottom_(pos))
                            info.set(Indices::velocityYIdx);     // no-slip (top/bottom)
                        constraints_.push_back(DirichletConstraintData{
                            std::move(info), ConstraintValues(0.0), localDof.dofIndex()});
                    }
                }
            }
            else
            {
                static const bool anchor = getParam<bool>("Problem.AnchorPressureAtBottom", true);
                if (anchor)
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        const auto& pos = scv.dofPosition();
                        if (pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + 1e-6
                            && pos[0] < this->gridGeometry().bBoxMin()[0] + 1e-6)
                        {
                            ConstraintInfo info;
                            info.set(Indices::pressureIdx);
                            constraints_.push_back(DirichletConstraintData{
                                std::move(info), ConstraintValues(0.0), scv.dofIndex()});
                        }
                    }
            }
        }
    }

    Scalar surfaceTension_, interfaceThickness_, mobility_, scaledSurfaceTension_, energyScale_;
    Scalar rho1_, rho2_, eta1_, eta2_, g_, radius_;
    SurfaceTensionForm surfaceTensionForm_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
