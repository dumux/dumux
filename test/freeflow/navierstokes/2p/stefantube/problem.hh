// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_TWOP_STEFAN_TUBE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_TWOP_STEFAN_TUBE_PROBLEM_HH

#include <algorithm>
#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \brief Slit-capillary Stefan-tube evaporation benchmark.
 *
 * Liquid is on the left, gas is on the right, and the opening at x=L fixes
 * the far-field vapor concentration. The meniscus is initialized as a
 * circular arc with prescribed contact angle. The analytical reference is the
 * cross-section-averaged tube law ell^2 = ell0^2 + 2Kt.
 */
template <class TypeTag, class BaseProblem>
class StefanTubeProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using DirichletValues = typename ParentType::DirichletValues;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Vertex = typename GridGeometry::GridView::template Codim<dim>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    StefanTubeProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension", 1e-3);
        interfaceThickness_ = getParam<Scalar>("Problem.InterfaceThickness", 5e-3);
        mobility_ = getParam<Scalar>("Problem.Mobility", 5e-6);

        const auto gamma = surfaceTension_ * 3.0/(2.0*std::sqrt(2.0));
        scaledSurfaceTension_ = gamma;
        surfaceTensionCoeff_ = gamma*interfaceThickness_;
        energyScale_ = gamma/interfaceThickness_;

        rhoLiquid_ = getParam<Scalar>("Problem.Density1", 1.0);
        rhoGas_ = getParam<Scalar>("Problem.Density2", 1e-2);
        etaLiquid_ = getParam<Scalar>("Problem.Viscosity1", 1.0);
        etaGas_ = getParam<Scalar>("Problem.Viscosity2", 1e-2);

        tubeMinX_ = this->gridGeometry().bBoxMin()[0];
        openingX_ = this->gridGeometry().bBoxMax()[0];
        if constexpr (dimWorld == 3)
        {
            // The 3D grid covers only one quarter of the tube cross-section: the
            // bBoxMin sides in y and z are the two symmetry (centerline) planes,
            // and only the bBoxMax sides are the true, physical tube wall. The
            // full tube height is therefore twice the domain's own y-extent, and
            // the true centerline sits at the domain's y=bBoxMin (not its middle).
            tubeHeight_ = 2.0*(this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]);
            tubeCenterY_ = this->gridGeometry().bBoxMin()[1];
            tubeCenterZ_ = this->gridGeometry().bBoxMin()[2];
        }
        else
        {
            tubeHeight_ = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
            tubeCenterY_ = 0.5*(this->gridGeometry().bBoxMax()[1] + this->gridGeometry().bBoxMin()[1]);
            tubeCenterZ_ = 0.0;
        }

        const auto defaultGasLength = openingX_ - getParam<Scalar>("Problem.InterfacePosition", 0.70);
        gasLength0_ = getParam<Scalar>("Problem.InitialGasLength", defaultGasLength);
        gasLength0_ = std::clamp(gasLength0_, Scalar(10.0)*interfaceThickness_, openingX_ - tubeMinX_);
        flatInitialInterface_ = getParam<bool>("Problem.FlatInitialInterface", false);
        // Corner-filament scenario (3D only): below the Concus-Finn threshold
        // (ContactAngle < 45 deg for this square corner) a corner filament can
        // wick along the entire available wall length; this seeds the IC as
        // already wicked to near the opening in the corner region (where the
        // bulk spherical cap is invalid, dy^2+dz^2 > arcRadius^2), instead of
        // flat-clamping it at the bulk interface position. See
        // corner-meniscus-notes.md.
        cornerFilamentIC_ = getParam<bool>("Problem.CornerFilamentInitialCondition", false);

        enableEvaporation_ = getParam<bool>("Problem.EnableEvaporation", true);
        pinMomentum_ = getParam<bool>("Problem.PinMomentum", false);
        enableCapillaryMomentumSource_ = getParam<bool>("Problem.EnableCapillaryMomentumSource", true);
        evaporationRateCoeff_ = getParam<Scalar>("Problem.EvaporationRateCoeff", 500.0);
        vaporSatConc_ = getParam<Scalar>("Problem.VaporSatConc", 0.20);
        vaporFarFieldConc_ = getParam<Scalar>("Problem.VaporFarFieldConc", 0.05);
        vaporDiffusivity_ = getParam<Scalar>("Problem.VaporDiffusivity", 1e-2);
        vaporInterfaceConc_ = vaporSatConc_*getParam<Scalar>("Problem.KelvinFactor", 1.0);

        phaseChangeFactor_ = getParam<Scalar>("Problem.PhaseChangeFactor", 2.0);

        const auto contactAngleDeg = getParam<Scalar>("Problem.ContactAngle", 60.0);
        contactAngle_ = contactAngleDeg*M_PI/180.0;
        const auto cosContactAngle = std::cos(contactAngle_);
        const auto absCos = std::abs(cosContactAngle);
        const auto defaultConcavity = cosContactAngle >= 0.0 ? -1.0 : 1.0;
        arcSign_ = getParam<Scalar>("Problem.MeniscusConcavity", defaultConcavity) >= 0.0 ? 1.0 : -1.0;
        if (absCos > 1e-8)
        {
            arcRadius_ = tubeHeight_/(2.0*absCos);
            arcMeanSqrt_ = meanArcSqrt_();
        }
        else
        {
            arcRadius_ = 1e30;
            arcMeanSqrt_ = arcRadius_;
            arcSign_ = 0.0;
        }

        tubeCoefficient_ = vaporDiffusivity_
            *std::max(Scalar(0), vaporInterfaceConc_ - vaporFarFieldConc_)
            /rhoLiquid_;
    }

    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolume& scv) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto& pos = scv.dofPosition();
            bool onSymmetryPlane = false;
            if constexpr (dimWorld == 3)
                for (int dir = 1; dir < dimWorld; ++dir)
                    onSymmetryPlane = onSymmetryPlane || isSymmetryPlane_(pos, dir);

            if (pinMomentum_)
                values.setAllDirichlet();
            else if (isRight_(pos))
                values.setAllNeumann();
            else if (isSolidWall_(pos))
                // The true tube wall wins even where it meets a symmetry plane
                // (e.g. y=bBoxMax, z=0): full no-slip in all directions, not
                // free-slip. Without this check first, onSymmetryPlane below is
                // triggered by the z=0 membership alone and wrongly gives this
                // edge free-slip/Neumann tangential momentum -- confirmed via a
                // large, sharply localized pressure spike (order 10, vs. a
                // domain-wide std of ~0.7) sitting exactly on this edge.
                values.setAllDirichlet();
            else if (onSymmetryPlane && !isLeft_(pos))
            {
                // Quarter-tube symmetry plane (not the physical tube wall): zero
                // normal velocity (no penetration), free-slip/natural (Neumann)
                // tangential momentum. Standard DuMux symmetry-axis pattern, see
                // e.g. test/freeflow/navierstokes/channel/pipe/momentum/problem.hh.
                for (int dir = 0; dir < dimWorld; ++dir)
                {
                    if (dir >= 1 && isSymmetryPlane_(pos, dir))
                        values.setDirichlet(Indices::velocity(dir));
                    else
                        values.setNeumann(Indices::momentumBalanceIdx(dir));
                }
            }
            else
                values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
            const auto& pos = scv.dofPosition();

            if (isLeft_(pos) || isRight_(pos))
                values.setDirichlet(Indices::phaseFieldIdx);

            if (isRight_(pos))
                values.setDirichlet(Indices::vaporIdx);

        }

        return values;
    }

    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        DirichletValues values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            if (isLeft_(globalPos))
                values[Indices::phaseFieldIdx] = 1.0;
            else if (isRight_(globalPos))
                values[Indices::phaseFieldIdx] = -1.0;

            if (isRight_(globalPos))
                values[Indices::vaporIdx] = vaporFarFieldConc_;
        }

        return values;
    }

    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        Sources source(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto& volVars = elemVolVars[scv];
            const auto phi = volVars.phaseField();

            source[Indices::chemicalPotentialEqIdx] =
                volVars.chemicalPotential() - energyScale_*phi*(phi*phi - 1.0);

            if (enableEvaporation_)
            {
                const auto mDot = evaporationMassSource(phi, volVars.vaporConcentration());
                source[Indices::phaseFieldEqIdx] -= phaseChangeFactor_*mDot/rhoLiquid_;
                source[Indices::vaporEqIdx] += mDot;
            }
        }
        else if (enableCapillaryMomentumSource_)
        {
            const auto ip = ipData(fvGeometry, scv);
            const auto grads = this->couplingManager().gradients(element, fvGeometry, ip);
            const auto values = this->couplingManager().values(element, fvGeometry, ip);

            const auto mu = values[CouplingManager::chemicalPotentialIdx];
            source = mu*grads[CouplingManager::phaseFieldIdx];
        }

        return source;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const typename FVElementGeometry::SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            if (isRight_(scvf.ipGlobal()))
                values[Indices::conti0EqIdx] =
                    this->faceVelocity(element, fvGeometry, scvf)*scvf.unitOuterNormal();

            if (isSolidWall_(scvf.ipGlobal()))
            {
                // Contact angle BC on the tube wall, measured through the liquid phase:
                //   gamma*eps*(grad phi . n) = f_w'(phi),
                //   f_w'(phi) = -(gamma cos(theta))*(3/4)*(1 - phi^2).
                // Here gamma is the scaled surface tension and phi=+1 denotes liquid.
                const auto phi = elemVolVars[scvf.insideScvIdx()].phaseField();
                values[Indices::chemicalPotentialEqIdx] =
                    -scaledSurfaceTension_*std::cos(contactAngle_)*(3.0/4.0)*(1.0 - phi*phi);
            }
        }

        return values;
    }

    InitialValues initialAtPos(const GlobalPosition& pos) const
    {
        InitialValues values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            values[Indices::pressureIdx] = 0.0;

            auto xi = interfacePositionAt(pos, 0.0);
            if constexpr (dimWorld == 3)
            {
                if (cornerFilamentIC_)
                {
                    const auto dy = pos[1] - tubeCenterY_;
                    const auto dz = pos[2] - tubeCenterZ_;
                    if (dy*dy + dz*dz > arcRadius_*arcRadius_)
                        xi = std::max(xi, openingX_ - Scalar(2.0)*interfaceThickness_);
                }
            }

            values[Indices::phaseFieldIdx] =
                std::tanh((xi - pos[0])/(std::sqrt(2.0)*interfaceThickness_));
            values[Indices::chemicalPotentialIdx] = 0.0;
            // Seed with the quasi-steady linear profile itself (the same one the tube-law
            // reference assumes at every instant) rather than an ad hoc plateau-plus-ramp guess.
            // A mismatched IC otherwise drives a large early-time transient in the evaporation
            // rate (observed: >3x overshoot at very early t) before it settles.
            values[Indices::vaporIdx] = analyticalVaporConcentration(pos, 0.0);
        }

        return values;
    }

    InitialValues initial(const SubControlVolume& scv) const
    { return initialAtPos(scv.dofPosition()); }

    InitialValues initial(const Vertex& vertex) const
    { return initialAtPos(vertex.geometry().corner(0)); }

    template<class Context>
    auto interfaceFlux(const Context& context) const
    {
        Dune::FieldVector<Scalar, dimWorld> zero(0.0);
        return zero;
    }

    // Used by the optional corner-wicking diagnostic in main.cc: forces the
    // evaporation source to zero so only Cahn-Hilliard/wetting relaxation acts.
    void setRelaxationPhase(bool value)
    { relaxationPhase_ = value; }

    Scalar mobility() const
    { return mobility_; }

    Scalar surfaceTension() const
    { return surfaceTensionCoeff_; }

    Scalar vaporDiffusivity() const
    { return vaporDiffusivity_; }

    Scalar mixtureViscosity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, Scalar(-1.0), Scalar(1.0));
        return 0.5*((1.0 + phi)*etaLiquid_ + (1.0 - phi)*etaGas_);
    }

    Scalar mixtureDensity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, Scalar(-1.0), Scalar(1.0));
        return 0.5*((1.0 + phi)*rhoLiquid_ + (1.0 - phi)*rhoGas_);
    }

    Scalar mixtureDensityDerivative() const
    { return 0.5*(rhoLiquid_ - rhoGas_); }

    Scalar liquidVolumeFraction(const Scalar phaseField) const
    { return 0.5*(1.0 + std::clamp(phaseField, Scalar(-1.0), Scalar(1.0))); }

    Scalar gasVolumeFraction(const Scalar phaseField) const
    { return 1.0 - liquidVolumeFraction(phaseField); }

    Scalar interfacialDelta(const Scalar phaseField) const
    {
        const auto g = std::max(Scalar(0), 1.0 - phaseField*phaseField);
        return g/(2.0*std::sqrt(2.0)*interfaceThickness_);
    }

    Scalar evaporationMassSource(const Scalar phaseField, const Scalar vaporConcentration) const
    {
        if (!enableEvaporation_ || relaxationPhase_)
            return 0.0;

        const auto drivingForce = std::max(Scalar(0), vaporInterfaceConc_ - vaporConcentration);
        return evaporationRateCoeff_*drivingForce*interfacialDelta(phaseField);
    }

    Scalar analyticalGasLength(const Scalar time) const
    { return std::sqrt(gasLength0_*gasLength0_ + 2.0*tubeCoefficient_*time); }

    Scalar analyticalInterfacePosition(const Scalar time) const
    { return openingX_ - analyticalGasLength(time); }

    // Meniscus surface: a circular arc in 2D (curved in y only, the exact 2D-slit
    // wetting solution), or a spherical cap in 3D (curved in y AND z jointly, using
    // the SAME radius arcRadius_ = tubeHeight_/(2|cos theta|) -- by symmetry a sphere
    // of this radius, centered on the tube axis, meets ALL FOUR walls of a square
    // tube at exactly the prescribed contact angle, away from the corner edges
    // themselves). A y-only profile would only satisfy the wetting condition
    // (gamma_eps grad(phi).n = ...) at the y-walls: its z-gradient is identically
    // zero, mismatching the nonzero z-wall condition for theta != 90 deg.
    Scalar interfacePositionAt(const GlobalPosition& pos, const Scalar time) const
    {
        if (flatInitialInterface_)
            return analyticalInterfacePosition(time);

        const auto dy = pos[1] - tubeCenterY_;
        if constexpr (dimWorld == 3)
        {
            const auto dz = pos[2] - tubeCenterZ_;
            const auto radicand = std::max(Scalar(0), arcRadius_*arcRadius_ - dy*dy - dz*dz);
            return analyticalInterfacePosition(time)
                + arcSign_*(std::sqrt(radicand) - arcMeanSqrt_);
        }
        else
        {
            const auto radicand = std::max(Scalar(0), arcRadius_*arcRadius_ - dy*dy);
            return analyticalInterfacePosition(time)
                + arcSign_*(std::sqrt(radicand) - arcMeanSqrt_);
        }
    }

    Scalar analyticalMassFlux(const Scalar time) const
    {
        const auto ell = analyticalGasLength(time);
        return ell > 0.0 ? rhoLiquid_*tubeCoefficient_/ell : 0.0;
    }

    Scalar analyticalVaporConcentration(const GlobalPosition& pos, const Scalar time) const
    {
        const auto xi = interfacePositionAt(pos, time);
        const auto ell = openingX_ - xi;
        if (ell <= 0.0)
            return vaporFarFieldConc_;

        const auto s = std::clamp((pos[0] - xi)/ell, Scalar(0), Scalar(1));
        return vaporInterfaceConc_ + (vaporFarFieldConc_ - vaporInterfaceConc_)*s;
    }

    Scalar openingPosition() const
    { return openingX_; }

    Scalar tubeHeight() const
    { return tubeHeight_; }

    Scalar liquidDensity() const
    { return rhoLiquid_; }

    Scalar vaporInterfaceConc() const
    { return vaporInterfaceConc_; }

    Scalar vaporFarFieldConc() const
    { return vaporFarFieldConc_; }

    Scalar tubeCoefficient() const
    { return tubeCoefficient_; }

private:
    bool isLeft_(const GlobalPosition& pos) const
    { return pos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isRight_(const GlobalPosition& pos) const
    { return pos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isBottom_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps_; }

    bool isTop_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_; }

    bool isSolidWall_(const GlobalPosition& pos) const
    {
        for (int dir = 1; dir < dimWorld; ++dir)
        {
            if constexpr (dimWorld == 3)
            {
                // bBoxMin sides are the quarter-tube symmetry planes, not walls.
                if (pos[dir] > this->gridGeometry().bBoxMax()[dir] - eps_)
                    return true;
            }
            else
            {
                if (pos[dir] < this->gridGeometry().bBoxMin()[dir] + eps_
                    || pos[dir] > this->gridGeometry().bBoxMax()[dir] - eps_)
                    return true;
            }
        }
        return false;
    }

    // Quarter-tube symmetry plane at the bBoxMin side of transverse direction dir
    // (3D only; the 2D domain is not symmetry-reduced and this is never true there).
    bool isSymmetryPlane_(const GlobalPosition& pos, int dir) const
    { return pos[dir] < this->gridGeometry().bBoxMin()[dir] + eps_; }

    // Cross-section area-average of the meniscus radicand: a 1D line-average of the
    // circular-arc profile in 2D, or a 2D area-average of the spherical-cap profile
    // over the full (y,z) cross-section in 3D (see interfacePositionAt).
    Scalar meanArcSqrt_() const
    {
        if constexpr (dimWorld == 3)
        {
            const int n = 256;
            Scalar sum = 0.0;
            for (int i = 0; i <= n; ++i)
            {
                const auto y = -0.5*tubeHeight_ + tubeHeight_*Scalar(i)/Scalar(n);
                const auto wy = (i == 0 || i == n) ? 0.5 : 1.0;
                for (int j = 0; j <= n; ++j)
                {
                    const auto z = -0.5*tubeHeight_ + tubeHeight_*Scalar(j)/Scalar(n);
                    const auto wz = (j == 0 || j == n) ? 0.5 : 1.0;
                    sum += wy*wz*std::sqrt(std::max(Scalar(0), arcRadius_*arcRadius_ - y*y - z*z));
                }
            }
            return sum/Scalar(n*n);
        }
        else
        {
            const int n = 512;
            Scalar sum = 0.0;
            for (int i = 0; i <= n; ++i)
            {
                const auto y = -0.5*tubeHeight_ + tubeHeight_*Scalar(i)/Scalar(n);
                const auto w = (i == 0 || i == n) ? 0.5 : 1.0;
                sum += w*std::sqrt(std::max(Scalar(0), arcRadius_*arcRadius_ - y*y));
            }
            return sum/Scalar(n);
        }
    }

    static constexpr Scalar eps_ = 1e-8;

    Scalar surfaceTension_;
    Scalar interfaceThickness_;
    Scalar mobility_;
    Scalar scaledSurfaceTension_;
    Scalar surfaceTensionCoeff_;
    Scalar energyScale_;

    Scalar rhoLiquid_;
    Scalar rhoGas_;
    Scalar etaLiquid_;
    Scalar etaGas_;

    Scalar tubeMinX_;
    Scalar openingX_;
    Scalar tubeHeight_;
    Scalar tubeCenterY_;
    Scalar tubeCenterZ_;
    Scalar gasLength0_;
    bool flatInitialInterface_;
    bool cornerFilamentIC_;
    Scalar contactAngle_;
    Scalar arcRadius_;
    Scalar arcMeanSqrt_;
    Scalar arcSign_;

    bool enableEvaporation_;
    bool relaxationPhase_ = false;
    bool pinMomentum_;
    bool enableCapillaryMomentumSource_;
    Scalar evaporationRateCoeff_;
    Scalar vaporSatConc_;
    Scalar vaporInterfaceConc_;
    Scalar vaporFarFieldConc_;
    Scalar vaporDiffusivity_;
    Scalar phaseChangeFactor_;
    Scalar tubeCoefficient_;
};

} // end namespace Dumux

#endif
