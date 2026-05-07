// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Momentum sub-problem for the gripper hydrogel actuator.
 *
 * Provides boundary conditions, initial conditions, material parameters,
 * and the optional rigid-sphere contact source term for the momentum
 * subdomain (displacement \f$\mathbf{u}\f$, PQ1Bubble discretisation).
 *
 * \par Boundary conditions
 * - Left end (\f$x = x_\mathrm{min}\f$): clamped, \f$\mathbf{u} = \mathbf{0}\f$ (Dirichlet).
 * - All other boundaries: traction-free (homogeneous Neumann).
 *
 * \par Rigid-sphere contact (optional)
 * If \c Problem.EnablePearlContact is set, a quadratic penalty force is added
 * as a source term whenever a DOF penetrates the sphere of radius \f$R_p\f$
 * centred at \f$\mathbf{c}\f$:
 * \f[
 *   \mathbf{f}_\mathrm{contact}
 *   = \varepsilon_p \bigl(R_p - \|\mathbf{x}_\mathrm{cur} - \mathbf{c}\|\bigr)^2
 *     \hat{\mathbf{n}}_\mathrm{out}
 * \f]
 * The quadratic law ensures a continuous Jacobian at first contact, improving
 * Newton convergence compared to a linear penalty.
 *
 * \par Dynamic time integration
 * Inertia \f$\rho\,\ddot{\mathbf{u}}\f$ and optional velocity damping
 * \f$\eta\,\dot{\mathbf{u}}\f$ are evaluated via a Newmark-beta scheme
 * injected through \c setNewmarkScheme().
 */
#ifndef DUMUX_GRIPPER_PROBLEM_MOMENTUM_HH
#define DUMUX_GRIPPER_PROBLEM_MOMENTUM_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/experimental/timestepping/newmarkbeta.hh>

namespace Dumux {

template<class TypeTag>
class GripperMomentumProblem : public FVProblemWithSpatialParams<TypeTag>
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
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;

public:
    GripperMomentumProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                           std::shared_ptr<CouplingManager> couplingManager,
                           const std::string& paramGroup = "Momentum")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        xMin_ = this->gridGeometry().bBoxMin()[0];
        xMax_ = this->gridGeometry().bBoxMax()[0];
        yMin_ = this->gridGeometry().bBoxMin()[1];
        yMax_ = this->gridGeometry().bBoxMax()[1];
        if constexpr (dim > 2)
        {
            zMin_ = this->gridGeometry().bBoxMin()[2];
            zMax_ = this->gridGeometry().bBoxMax()[2];
        }

        enablePearlContact_ = getParam<bool>("Problem.EnablePearlContact", false);
        if (enablePearlContact_)
        {
            pearlRadius_   = getParam<Scalar>("Problem.PearlRadius");
            pearlCenter_[0] = getParam<Scalar>("Problem.PearlCenterX");
            pearlCenter_[1] = getParam<Scalar>("Problem.PearlCenterY");
            if constexpr (dim > 2)
                pearlCenter_[2] = getParam<Scalar>("Problem.PearlCenterZ", 0.5*(zMin_ + zMax_));
            contactPenalty_ = getParam<Scalar>("Problem.ContactPenalty");
        }

        solidDensity_ = getParam<Scalar>("SpatialParams.SolidDensity", Scalar(1000.0));
        viscousDamping_ = getParam<Scalar>("Momentum.ViscousDamping", Scalar(0.0));
    }

    //! Return the coupling manager.
    const CouplingManager& couplingManager() const { return *couplingManager_; }

    //! Set the Newmark-beta scheme used for dynamic time integration.
    void setNewmarkScheme(std::shared_ptr<const Experimental::NewmarkBeta<Scalar, SolutionVector>> newmark)
    { newmark_ = std::move(newmark); }

    //! Effective solid density used in inertia term.
    Scalar solidDensity(const Element& /*element*/, const SubControlVolume& /*scv*/) const
    { return solidDensity_; }

    //! Velocity-proportional damping coefficient for transient momentum equation.
    Scalar viscousDamping() const
    { return viscousDamping_; }

    //! Compute acceleration from Newmark-beta at current Newton iterate displacement.
    auto acceleration(const Element& /*element*/,
                      const SubControlVolume& scv,
                      const Scalar dt,
                      const Dune::FieldVector<Scalar, dim>& displacement) const
    {
        return newmark_->acceleration(scv.dofIndex(), dt, displacement);
    }

    /*!
     * \brief Boundary types.
     *
     * The left end (x ≈ x_min) is clamped.  Everything else is Neumann.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (onLeftEnd_(globalPos))
            values.setAllDirichlet();
        return values;
    }

    //! Dirichlet values: zero displacement (clamped).
    PrimaryVariables dirichletAtPos(const GlobalPosition& /*globalPos*/) const
    { return PrimaryVariables(0.0); }

    //! Neumann fluxes: traction-free (zero).
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& /*element*/,
                        const FVElementGeometry& /*fvGeometry*/,
                        const ElementVolumeVariables& /*elemVolVars*/,
                        const ElementFluxVariablesCache& /*elemFluxVarsCache*/,
                        const SubControlVolumeFace& /*scvf*/) const
    { return NumEqVector(0.0); }

    //! Initial condition: zero displacement.
    PrimaryVariables initialAtPos(const GlobalPosition& /*globalPos*/) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Optional pearl contact source term.
     *
     * Adds a penalty force on any DOF whose current position has penetrated
     * inside the analytical sphere of radius \f$R_p\f$:
     * \f[
    *   \mathbf{f}_i = \epsilon_p\, (-g)\, \hat{n}_\text{outward},
     *   \quad g = \|\mathbf{x}_\text{cur} - \mathbf{c}\| - R_p < 0.
     * \f]
    * A linear penalty response is stiffer near contact onset than a quadratic
    * law and helps prevent visible tunneling for coarse dynamic steps.
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector src(0.0);
        if (!enablePearlContact_) return src;

        // Current deformed position of this DOF.
        GlobalPosition curPos = scv.dofPosition();
        for (int dir = 0; dir < dim; ++dir)
            curPos[dir] += elemVolVars[scv].displacement(dir);

        const GlobalPosition diff = curPos - pearlCenter_;
        const Scalar dist = diff.two_norm();
        const Scalar gap = dist - pearlRadius_;

        if (gap < 0.0)
        {
            // Linear penalty: f = ε·(-gap). Stronger contact onset response.
            const GlobalPosition nOut = diff / std::max(dist, Scalar(1e-12));
            const Scalar f = contactPenalty_ * (-gap);
            for (int dir = 0; dir < dim; ++dir)
                src[Indices::momentum(dir)] = f * nOut[dir];
        }
        return src;
    }

    //! The problem name (used for VTK output).
    std::string name() const { return "gripper_momentum"; }

private:
    bool onLeftEnd_(const GlobalPosition& p) const
    { return p[0] < xMin_ + eps_ * (xMax_ - xMin_); }

    static constexpr Scalar eps_ = 1e-6;
    Scalar xMin_, xMax_, yMin_, yMax_, zMin_ = 0.0, zMax_ = 0.0;

    bool enablePearlContact_ = false;
    Scalar pearlRadius_ = 0.0;
    GlobalPosition pearlCenter_{0.0};
    Scalar contactPenalty_ = 0.0;
    Scalar solidDensity_ = 1000.0;
    Scalar viscousDamping_ = 0.0;

    std::shared_ptr<const Experimental::NewmarkBeta<Scalar, SolutionVector>> newmark_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif // DUMUX_GRIPPER_PROBLEM_MOMENTUM_HH
