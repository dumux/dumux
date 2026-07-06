// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawTests
 * \brief Saffman-Taylor fingering test for the hybrid CVFE Hele-Shaw
 *        Darcy-Cahn-Hilliard model.
 *
 * Domain: [0, L] x [0, H], flow from left (high pressure) to right (low pressure).
 * Less viscous fluid (phi=+1) displaces more viscous fluid (phi=-1).
 * The interface is initially perturbed to trigger fingering.
 *
 * Pressure is essential (Dirichlet) at the inlet/outlet, imposed via
 * `constraints()`; every other boundary condition (the no-flow walls, and the
 * advective outflow of the phase field at the inlet/outlet) is weakly
 * imposed via `boundaryFlux()`.
 */
#ifndef DUMUX_TEST_PHASEFIELD_HELESHAW_2P_CVFE_PROBLEM_HH
#define DUMUX_TEST_PHASEFIELD_HELESHAW_2P_CVFE_PROBLEM_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/concepts/ipdata_.hh>

#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>

#include <dumux/phasefield/heleshaw/2p/cvfe/flux.hh>
#include <dumux/phasefield/freeenergy/symmetricdoublewell.hh>

namespace Dumux {

/*!
 * \ingroup HeleShawTests
 * \brief Test problem for the hybrid CVFE Hele-Shaw model on PQ2 elements.
 */
template<class TypeTag>
class HeleShawTwoPCVFETestProblem : public Dumux::Experimental::Problem<TypeTag>
{
    using ParentType = Dumux::Experimental::Problem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;

    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables, GridIndexType>;

    using GravityVector = Dune::FieldVector<Scalar, GridView::dimensionworld>;

public:
    HeleShawTwoPCVFETestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , freeEnergy_(computeEnergyScale_(getParam<Scalar>("Problem.SurfaceTension"),
                                      getParam<Scalar>("Problem.InterfaceThickness")))
    {
        pInlet_ = getParam<Scalar>("Problem.PressureInlet");
        pOutlet_ = getParam<Scalar>("Problem.PressureOutlet", 0.0);

        rho1_ = getParam<Scalar>("Problem.DensityInvading", 1000.0);
        rho2_ = getParam<Scalar>("Problem.DensityReceding",  1000.0);
        eta1_ = getParam<Scalar>("Problem.ViscosityInvading", 1.0);
        eta2_ = getParam<Scalar>("Problem.ViscosityReceding", 10.0);

        permeability_ = getParam<Scalar>("Problem.Permeability");

        const auto sigma = getParam<Scalar>("Problem.SurfaceTension");
        epsilon_ = getParam<Scalar>("Problem.InterfaceThickness");
        chMobility_ = getParam<Scalar>("Problem.CHMobility", epsilon_*epsilon_);

        // gamma = 3*sigma/(2*sqrt(2)), surfaceTensionCoeff = gamma*epsilon
        const auto gamma = sigma*3.0/(2.0*std::sqrt(2.0));
        surfaceTensionCoeff_ = gamma*epsilon_;

        x0_ = getParam<Scalar>("Problem.InterfacePosition");
        delta_ = getParam<Scalar>("Problem.InterfacePerturbation", 0.0);
        nModes_ = getParam<int>("Problem.PerturbationModes", 1);
        baseMode_ = getParam<int>("Problem.PerturbationBaseMode", 1);

        const auto g = getParam<Scalar>("Problem.Gravity", 0.0);
        gravity_ = {0.0, -g};

        appendDirichletConstraints_();
    }

    // ------- boundary conditions -------

    /*!
     * \brief All equations are flux (weakly imposed) boundaries by default;
     *        the pressure Dirichlet condition at the inlet/outlet is imposed
     *        on top via `constraints()`, which fully overrides the residual
     *        and Jacobian row of the pressure equation at those local dofs.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    /*!
     * \brief Essential (Dirichlet) pressure constraints at the inlet/outlet.
     */
    const auto& constraints() const
    { return constraints_; }

    /*!
     * \brief Weakly imposed (Neumann) boundary fluxes.
     *
     * Top/bottom walls: no-flow for all equations. Inlet/outlet: the
     * advective outflow of phi is the key flux (without it, the boundary
     * dof has no outflux and phi creeps toward the wrong value); mu keeps
     * its natural (zero-diffusion) condition; the pressure component is
     * never used since pressure is Dirichlet there via `constraints()`.
     */
    template<class ElementVariables, class FaceIpData>
    NumEqVector boundaryFlux(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const FaceIpData& faceIpData) const
    {
        NumEqVector flux(0.0);
        const auto& pos = faceIpData.global();

        if (!isInlet_(pos) && !isOutlet_(pos))
            return flux;

        const HeleShawTwoPCVFEFluxFunctionContext context(fvGeometry, elemVars, cache(elemVars, faceIpData));

        const auto phi = context.phaseField();
        const auto mu = context.chemicalPotential();
        const auto etaMix = mixtureViscosity(phi);
        const auto rhoMix = mixtureDensity(phi);
        const auto lambda = permeability_/etaMix;
        const auto& n = faceIpData.unitOuterNormal();

        // Normal Darcy velocity (per unit area; assembly multiplies by area)
        const auto vn = -lambda*(context.gradPressure()*n - mu*(context.gradPhaseField()*n) - rhoMix*(gravity_*n));

        // Phase field: upwind from inside for outflow, prescribed BC value for inflow.
        // Inlet carries phi=+1 (invading fluid), outlet carries whatever is inside.
        const Scalar phiBC = isInlet_(pos) ? Scalar(1.0) : Scalar(-1.0);
        const auto phiUpwind = (vn < 0.0) ? phiBC : phi;

        flux[Indices::phaseFieldEqIdx] = phiUpwind*vn;
        return flux;
    }

    // ------- initial conditions -------

    PrimaryVariables initialAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(0.0);

        // Sinusoidally perturbed interface position x_int(y)
        const auto H = this->gridDiscretization().bBoxMax()[1];
        const auto L = this->gridDiscretization().bBoxMax()[0];
        auto xInt = x0_;
        for (int k = baseMode_; k < baseMode_ + nModes_; ++k)
            xInt += delta_*std::cos(k*2.0*M_PI*pos[1]/H);

        // Smooth tanh profile: phi=+1 left of interface (invading), phi=-1 right (receding)
        values[Indices::phaseFieldIdx] = std::tanh((xInt - pos[0])/(std::sqrt(2.0)*epsilon_));
        values[Indices::chemPotIdx] = 0.0;

        // Linear pressure from inlet to outlet as initial guess
        values[Indices::pressureIdx] = pInlet_ + (pOutlet_ - pInlet_)*(pos[0]/L);

        return values;
    }

    // ------- source term -------

    /*!
     * \brief Adds the derivative of the double-well free energy to the
     *        chemical potential equation. The model-intrinsic identity
     *        contribution (+mu) is already added by the local residual.
     */
    template<class ElementVariables, class IpData>
    NumEqVector source(const FVElementGeometry& fvGeometry,
                       const ElementVariables& elemVars,
                       const IpData& ipData) const
    {
        NumEqVector values(0.0);

        Scalar phi = 0.0;
        if constexpr (Dumux::Concept::LocalDofIpData<IpData>)
            phi = elemVars[ipData].phaseField();
        else
            phi = HeleShawTwoPCVFEFluxFunctionContext(fvGeometry, elemVars, cache(elemVars, ipData)).phaseField();

        values[Indices::chemPotEqIdx] = -freeEnergy_.derivative(phi);
        return values;
    }

    // ------- fluid properties (called from Variables::update and LocalResidual) -------

    Scalar mixtureDensity(Scalar phi) const
    {
        phi = std::clamp(phi, Scalar(-1.0), Scalar(1.0));
        return 0.5*((1.0 + phi)*rho1_ + (1.0 - phi)*rho2_);
    }

    Scalar mixtureViscosity(Scalar phi) const
    {
        phi = std::clamp(phi, Scalar(-1.0), Scalar(1.0));
        return 0.5*((1.0 + phi)*eta1_ + (1.0 - phi)*eta2_);
    }

    Scalar permeability() const { return permeability_; }
    Scalar chMobility() const { return chMobility_; }
    Scalar surfaceTensionCoeff() const { return surfaceTensionCoeff_; }
    GravityVector gravity() const { return gravity_; }

private:
    bool isInlet_(const GlobalPosition& pos) const
    { return pos[0] < this->gridDiscretization().bBoxMin()[0] + 1e-8; }

    bool isOutlet_(const GlobalPosition& pos) const
    { return pos[0] > this->gridDiscretization().bBoxMax()[0] - 1e-8; }

    static Scalar computeEnergyScale_(Scalar sigma, Scalar epsilon)
    {
        const auto gamma = sigma*3.0/(2.0*std::sqrt(2.0));
        return gamma/epsilon;
    }

    void appendDirichletConstraints_()
    {
        auto fvGeometry = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            fvGeometry.bind(element);

            for (const auto& boundaryFace : boundaryFaces(fvGeometry))
            {
                const auto& pos = boundaryFace.center();
                if (isInlet_(pos) || isOutlet_(pos))
                {
                    for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                    {
                        const auto& globalPos = ipData(fvGeometry, localDof).global();

                        ConstraintInfo info;
                        info.set(Indices::pressureIdx);

                        PrimaryVariables dirichletValues(0.0);
                        dirichletValues[Indices::pressureIdx] = isInlet_(globalPos) ? pInlet_ : pOutlet_;

                        constraints_.push_back(DirichletConstraintData{std::move(info), std::move(dirichletValues), localDof.dofIndex()});
                    }
                }
            }
        }
    }

    Scalar pInlet_, pOutlet_;
    Scalar rho1_, rho2_, eta1_, eta2_;
    Scalar permeability_;
    Scalar epsilon_, chMobility_, surfaceTensionCoeff_;
    Scalar x0_, delta_;
    int nModes_, baseMode_;
    GravityVector gravity_;
    Dumux::FreeEnergy::SymmetricDoubleWell<Scalar> freeEnergy_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
