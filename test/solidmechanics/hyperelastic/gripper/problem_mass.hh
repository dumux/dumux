// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Mass sub-problem for the gripper hydrogel actuator.
 *
 * Provides boundary conditions and initial conditions for the mass subdomain,
 * which carries two primary variables per Box/P1 DOF:
 * - \f$p_\mathrm{fluid}\f$ — fluid pore pressure (Darcy flow).
 * - \f$p_s\f$ — solid bulk pressure (mixed u-p constraint, no spatial BC required).
 *
 * \par Boundary conditions for \f$p_\mathrm{fluid}\f$
 * - Active top surface (\f$y = y_\mathrm{max}\f$, \f$x > x_\mathrm{active}\f$)
 *   once \f$t \geq t_\mathrm{swell}\f$:
 *   Dirichlet \f$p_\mathrm{fluid} = \Pi(\phi_0)\,r(t)\f$,
 *   simulating immersion in solvent with a smooth load ramp
 *   \f$r(t) = \min\!\left((t-t_\mathrm{swell})/t_\mathrm{ramp},\,1\right)\f$.
 * - Passive anchor region (\f$x < x_\mathrm{active}\f$) top surface and all
 *   other boundaries: no-flux (Neumann = 0).
 *
 * The osmotic driving pressure \f$\Pi\f$ follows Flory-Huggins theory
 * (see `GripperSpatialParams::osmoticPressure`), evaluated at \f$J=1\f$.
 *
 * \par Boundary conditions for \f$p_s\f$
 * No explicit BC is required anywhere: \f$p_s\f$ is determined pointwise
 * by the local bulk pressure constraint \f$\partial U/\partial J = p_s/\lambda\f$
 * enforced in `GripperMassLocalResidual::computeSource`.
 */
#ifndef DUMUX_GRIPPER_PROBLEM_MASS_HH
#define DUMUX_GRIPPER_PROBLEM_MASS_HH

#include <algorithm>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

template<class TypeTag>
class GripperMassProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    GripperMassProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                           std::shared_ptr<CouplingManager> couplingManager,
                           const std::string& paramGroup = "Mass")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        yMin_ = this->gridGeometry().bBoxMin()[1];
        yMax_ = this->gridGeometry().bBoxMax()[1];
        xMin_ = this->gridGeometry().bBoxMin()[0];
        xMax_ = this->gridGeometry().bBoxMax()[0];

        swellingActiveTime_ = getParam<Scalar>("Problem.SwellingActiveTime", Scalar(0.0));
        swellingRampTime_   = getParam<Scalar>("Problem.SwellingRampTime",   Scalar(1.0));
        activeRegionStart_  = getParam<Scalar>("Problem.ActiveRegionStart", Scalar(0.2));

        // Osmotic driving pressure (Flory-Huggins at phi = phi0).
        // Evaluated at J=1 for the BC value (reference state for initial immersion).
        osmoticPressureBC_ = getParam<Scalar>(
            "SpatialParams.OsmoticPressureBC",
            this->spatialParams().osmoticPressure(Scalar(1.0)));
    }

    //! Return the coupling manager.
    const CouplingManager& couplingManager() const { return *couplingManager_; }

    //! Record the current simulation time (called from main.cc before each time step).
    void setTime(Scalar t) { time_ = t; }

    //! Whether the swelling boundary condition is currently active.
    bool swellingActive() const { return time_ >= swellingActiveTime_; }

    //! Ramp factor in [0,1]: linearly increases from 0 at swelling onset to 1 after ramp time.
    Scalar rampFactor() const
    {
        if (!swellingActive()) return Scalar(0.0);
        return std::min((time_ - swellingActiveTime_) / swellingRampTime_, Scalar(1.0));
    }

    /*!
     * \brief Boundary types.
     *
     * Top surface (y ≈ y_max) when swelling is active → Dirichlet pressure.
     * Everywhere else → Neumann (no-flux).
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        // Osmotic BC covers the active region plus transition zone (spatialRamp_ handles the
        // value tapering to zero at the start; the passive anchor stays Neumann).
        // p_s has no explicit BC — determined pointwise by the local constraint.
        if (onTopSurface_(globalPos) && onActiveRegion_(globalPos) && swellingActive())
            values.setDirichlet(0); // index 0 = pressureIdx
        return values;
    }

    /*!
     * \brief Dirichlet value on the top surface.
     *
     * Sets pore pressure to the osmotic driving value \f$p = -\Pi(\phi_0)\f$,
     * so that the chemical potential at the surface equals that of the pure solvent.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        if (onTopSurface_(globalPos) && swellingActive())
            values[0] = osmoticPressureBC_ * rampFactor();
        return values;
    }

    //! Neumann fluxes: no-flux everywhere (sealed, or active surface handled by Dirichlet).
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& /*element*/,
                        const FVElementGeometry& /*fvGeometry*/,
                        const ElementVolumeVariables& /*elemVolVars*/,
                        const ElementFluxVariablesCache& /*elemFluxVarsCache*/,
                        const SubControlVolumeFace& /*scvf*/) const
    { return NumEqVector(0.0); }

    //! Initial condition: zero pore pressure (equilibrated with atmosphere).
    PrimaryVariables initialAtPos(const GlobalPosition& /*globalPos*/) const
    { return PrimaryVariables(0.0); }

    //! The problem name.
    std::string name() const { return "gripper_pressure"; }

private:
    bool onTopSurface_(const GlobalPosition& p) const
    { return p[1] > yMax_ - eps_ * (yMax_ - yMin_); }

    bool onActiveRegion_(const GlobalPosition& p) const
    { return p[0] - xMin_ > activeRegionStart_ * (xMax_ - xMin_); }

    static constexpr Scalar eps_ = 1e-6;
    Scalar yMin_, yMax_, xMin_, xMax_;
    Scalar swellingActiveTime_ = 0.0;
    Scalar swellingRampTime_   = 1.0;
    Scalar activeRegionStart_  = 0.2;
    Scalar osmoticPressureBC_  = 0.0;
    Scalar time_ = 0.0;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif // DUMUX_GRIPPER_PROBLEM_MASS_HH
