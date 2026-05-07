// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Hydrogel actuator (gripper) model — two coupled subdomains on the same grid.
 *
 * \par Subdomain 0 — Momentum (displacement \f$\mathbf{u}\f$, PQ1Bubble)
 * Dynamic momentum balance for a compressible Neo-Hookean porous solid:
 * \f[
 *   \rho\,\ddot{\mathbf{u}} + \eta\,\dot{\mathbf{u}}
 *   - \nabla_X \cdot \mathbf{P} = \mathbf{f}_\mathrm{contact}
 * \f]
 * with the mixed u-p first Piola-Kirchhoff stress
 * \f[
 *   \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T})
 *              + J(p_s - p_\mathrm{fluid})\,\mathbf{F}^{-T}
 * \f]
 * where \f$\mathbf{F} = \mathbf{I} + \nabla_X\mathbf{u}\f$, \f$J = \det\mathbf{F}\f$.
 * The full strain energy is
 * \f$\psi = \tfrac{\lambda}{2}(\ln J)^2 + \tfrac{\mu}{2}(I_1-3-2\ln J)\f$.
 *
 * \par Subdomain 1 — Fluid/solid pressures (\f$p_\mathrm{fluid}\f$, \f$p_s\f$, Box/P1)
 * Two equations sharing P1 DOFs (MINI element pair with PQ1Bubble displacement):
 *
 * <b>Fluid mass balance:</b>
 * \f[
 *   \frac{\partial J}{\partial t}
 *   - \nabla_X \cdot \left(\frac{J K(J)}{\mu_f}
 *       \mathbf{F}^{-1}\mathbf{F}^{-T} \nabla_X p_\mathrm{fluid}\right) = 0
 * \f]
 * with permeability \f$K(J) = K_0 \exp(M(J-1)/J)\f$.
 *
 * <b>Solid bulk pressure constraint</b> (local, no spatial flux):
 * \f[
 *   \frac{\partial U}{\partial J}(J) - \frac{p_s}{\lambda} = 0,
 *   \qquad U(J) = \tfrac{1}{2}(\ln J)^2
 * \f]
 * This replaces the explicit \f$\lambda\ln J\f$ term in the stress; \f$\lambda\f$
 * no longer appears in \f$\mathbf{P}\f$ directly.
 *
 * \par Coupling
 * Both subdomains share the same grid. `GripperCouplingManager` provides:
 * - \f$p_\mathrm{fluid}\f$ and \f$p_s\f$ at SCVF quadrature points to the momentum residual.
 * - \f$\mathbf{F}\f$ and \f$J^\mathrm{prev}\f$ to the mass residual.
 */
#ifndef DUMUX_GRIPPER_MODEL_HH
#define DUMUX_GRIPPER_MODEL_HH

#include <cmath>
#include <memory>
#include <vector>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/multidomain/couplingmanager.hh>

#include "couplingmanager.hh"

namespace Dumux {

// ============================================================
// 1. Momentum subdomain (displacement, PQ1Bubble)
// ============================================================

/*!
 * \brief Equation and primary-variable indices for the momentum subdomain.
 * \tparam dim Spatial dimension.
 *
 * Displacement components \f$u_i\f$ occupy indices 0 to dim-1.
 */
template<int dim>
struct GripperMomentumIndices
{
    //! Index of the i-th displacement primary variable.
    static constexpr int displacement(int i) { return i; }
    //! Index of the i-th momentum balance equation.
    static constexpr int momentum(int i) { return i; }
};

/*!
 * \brief Model traits for the gripper momentum subdomain.
 */
template<int dim>
struct GripperMomentumModelTraits
{
    using Indices = GripperMomentumIndices<dim>;
    static constexpr int numEq() { return dim; }
};

//! VolumeVariables traits for momentum.
template<class PV, class MT>
struct GripperMomentumVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

/*!
 * \brief Volume variables for the momentum subdomain.
 * Stores the displacement vector.
 */
template<class Traits>
class GripperMomentumVolumeVariables : public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;
    //! Return the i-th displacement component.
    Scalar displacement(int i) const { return this->priVar(Indices::displacement(i)); }
};

/*!
 * \brief Local residual for the gripper momentum subdomain.
 *
 * Implements the quasi-static momentum balance:
 * \f[ -\nabla_X \cdot \mathbf{P} = 0,\quad
 *     \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T})
 *               + J(p_s - p_\text{fluid})\,\mathbf{F}^{-T} \f]
 *
 * Both \f$p_s\f$ (solid bulk pressure) and \f$p_\text{fluid}\f$ (pore pressure) live in
 * the mass subdomain (Box/P1) and are obtained via the coupling manager.
 */
template<class TypeTag>
class GripperMomentumLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::evalStorage;

    /*!
     * \brief Evaluate inertial storage contribution \f$\rho\,\ddot{\mathbf{u}}\f$.
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        Dune::FieldVector<Scalar, dimWorld> d(0.0);
        Dune::FieldVector<Scalar, dimWorld> dPrev(0.0);
        for (int dir = 0; dir < dimWorld; ++dir)
        {
            d[dir] = curVolVars.displacement(dir);
            dPrev[dir] = prevVolVars.displacement(dir);
        }

        const Scalar dt = this->timeLoop().timeStepSize();
        const auto a = problem.acceleration(element, scv, dt, d);
        const Scalar rho = problem.solidDensity(element, scv);
        const Scalar eta = problem.viscousDamping();
        const Scalar volume = curVolVars.extrusionFactor()*Extrusion::volume(fvGeometry, scv);

        for (int dir = 0; dir < dimWorld; ++dir)
        {
            residual[scv.localDofIndex()][Indices::momentum(dir)] += rho * a[dir] * volume;
            if (eta > 0.0)
                residual[scv.localDofIndex()][Indices::momentum(dir)] += eta * (d[dir] - dPrev[dir]) / dt * volume;
        }
    }

    /*!
     * \brief Momentum flux \f$-\mathbf{P}\,\hat{N}\,\mathrm{d}A\f$.
     *
     * The bulk pressure constraint has no spatial flux (it is a local algebraic equation
     * handled entirely in computeSource).
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux(0.0);
        const auto& n = scvf.unitOuterNormal();

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& fluxVarCache = elemFluxVarsCache[qpData.ipData()];
            const auto& ipGlobal = fluxVarCache.ipGlobal();

            // Deformation gradient F = ∇_X u + I at this integration point.
            Tensor F(0.0);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& gradN = fluxVarCache.gradN(scv.localDofIndex());
                for (int dir = 0; dir < dimWorld; ++dir)
                    F[dir].axpy(elemVolVars[scv.localDofIndex()].displacement(dir), gradN);
            }
            for (int dir = 0; dir < dimWorld; ++dir)
                F[dir][dir] += 1.0;

            // Both pressures live in the mass subdomain (Box/P1); evaluate via coupling.
            const Scalar p_solid = problem.couplingManager().solidBulkPressureAtPoint(
                fvGeometry, ipGlobal);
            const Scalar p_fluid = problem.couplingManager().pressureAtPoint(
                fvGeometry, ipGlobal);

            const Tensor P = piolaStress_(problem, element, F, p_solid, p_fluid);
            for (int i = 0; i < dimWorld; ++i)
                flux[Indices::momentum(i)] -= (P[i] * n) * qpData.weight();
        }
        return flux;
    }

private:
    /*!
     * \brief First Piola-Kirchhoff stress for the u-p mixed Neo-Hookean formulation.
     *
     * \f[
     *   \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T})
     *              + J(p_s - p_\text{fluid})\,\mathbf{F}^{-T}
     *            = \mu\mathbf{F} + \bigl(J(p_s - p_\text{fluid}) - \mu\bigr)\mathbf{F}^{-T}
     * \f]
     *
     * \f$\lambda\f$ does not appear here; it enters only through the bulk pressure
     * constraint \f$\partial U/\partial J = p_s/\kappa\f$ with \f$\kappa = \lambda\f$.
     */
    static Tensor piolaStress_(const Problem& problem,
                               const Element& element,
                               const Tensor& F,
                               const Scalar p_solid,
                               const Scalar p_fluid)
    {
        const Scalar J  = F.determinant();
        const Scalar mu = problem.spatialParams().shearModulus(element);

        Tensor Finv = F; Finv.invert();
        const Tensor FinvT = transpose(Finv);

        Tensor P(0.0);
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                P[i][j] = mu * F[i][j] + (J * (p_solid - p_fluid) - mu) * FinvT[i][j];
        return P;
    }
};

// ============================================================
// 2. Pressure subdomain (pore pressure, Box/P1)
// ============================================================

/*!
 * \brief Equation and primary-variable indices for the mass subdomain.
 *
 * Two primary variables per Box DOF:
 *  - index 0: fluid pore pressure \f$p_\text{fluid}\f$
 *  - index 1: solid bulk pressure \f$p_s\f$ (mixed u-p formulation)
 */
struct GripperMassIndices
{
    static constexpr int pressureIdx         = 0; //!< Fluid pore pressure.
    static constexpr int solidBulkPressureIdx = 1; //!< Solid bulk pressure p_s.
    static constexpr int massEqIdx           = 0; //!< Fluid mass balance equation.
    static constexpr int bulkPressureEqIdx    = 1; //!< Bulk pressure constraint equation.
};

//! Model traits for the gripper mass subdomain.
struct GripperMassModelTraits
{
    using Indices = GripperMassIndices;
    static constexpr int numEq() { return 2; } // p_fluid + p_solid
};

//! VolumeVariables traits for pressure.
template<class PV, class MT>
struct GripperMassVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

/*!
 * \brief Volume variables for the mass subdomain.
 * Stores fluid pore pressure and solid bulk pressure.
 */
template<class Traits>
class GripperMassVolumeVariables : public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;
    //! Return the fluid pore pressure \f$p_\text{fluid}\f$.
    Scalar pressure() const { return this->priVar(Indices::pressureIdx); }
    //! Return the solid bulk pressure \f$p_s\f$.
    Scalar solidBulkPressure() const { return this->priVar(Indices::solidBulkPressureIdx); }
};

/*!
 * \brief Local residual for the gripper mass subdomain.
 *
 * Two equations:
 *
 * **Fluid mass balance** (equation 0, primary variable \f$p_\text{fluid}\f$):
 * \f[ \frac{\partial J}{\partial t}
 *     - \nabla_X \cdot \!\left(\frac{J K(J)}{\mu_f}
 *         \mathbf{F}^{-1}\mathbf{F}^{-T} \nabla_X p_\text{fluid}\right) = 0 \f]
 *
 * **Solid bulk pressure constraint** (equation 1, primary variable \f$p_s\f$,
 * local — no spatial flux):
 * \f[ \frac{\partial U(J)}{\partial J} - \frac{p_s}{\kappa} = 0,\quad
 *     U(J) = \tfrac{1}{2}(\ln J)^2,\quad \kappa = \lambda \f]
 *
 * This P1 discretisation of \f$p_s\f$ paired with the PQ1Bubble displacement
 * satisfies the inf-sup condition (MINI element).
 */
template<class TypeTag>
class GripperMassLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;
    using ParentType::evalStorage;

    /*!
     * \brief Fluid mass storage: \f$(J^{n+1} - J^n) / \Delta t \cdot V_\text{scv}\f$.
     *
     * \f$J^{n+1}\f$ is obtained from the current momentum solution via the coupling manager.
     * \f$J^n\f$ is stored in the coupling manager and updated from main.cc after each
     * converged time step.
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        // J at the current Newton iterate from the coupled momentum domain.
        const Tensor F_cur = problem.couplingManager().deformationGradientAtPoint(fvGeometry, scv.dofPosition());
        const Scalar J_cur = F_cur.determinant();

        // J at the previous converged time step (stored in coupling manager).
        const Scalar J_prev = problem.couplingManager().prevJacobian(scv.dofIndex());

        const Scalar volume = curElemVolVars[scv].extrusionFactor()
                              * Extrusion::volume(fvGeometry, scv);
        const Scalar dt = this->timeLoop().timeStepSize();

        residual[scv.localDofIndex()][GripperMassIndices::massEqIdx]
            += (J_cur - J_prev) / dt * volume;
    }

    /*!
     * \brief Bulk pressure constraint: \f$\partial U(J)/\partial J - p_s/\kappa = 0\f$.
     *
     * Purely local — no spatial flux. \f$J = \det(\mathbf{F})\f$ is obtained from the
     * momentum domain via the coupling manager.
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        auto source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);

        const Tensor F = problem.couplingManager().deformationGradientAtPoint(
            fvGeometry, scv.dofPosition());
        const Scalar J    = F.determinant();
        const Scalar lnJ  = std::log(std::max(J, Scalar(1e-10)));
        const Scalar dUdJ = lnJ / std::max(J, Scalar(1e-10)); // U(J) = (lnJ)²/2
        const Scalar kappa = problem.spatialParams().firstLame(element); // κ = λ
        const Scalar p_s  = elemVolVars[scv].solidBulkPressure();

        source[GripperMassIndices::bulkPressureEqIdx] += dUdJ - p_s / kappa;
        return source;
    }

    /*!
     * \brief Darcy flux in the reference configuration.
     *
     * \f[
     *   q_\text{ref} = -\frac{J K(J)}{\mu_f}
     *                   \mathbf{F}^{-1}\mathbf{F}^{-T} \nabla_X p \cdot \hat{N}\,\mathrm{d}A
     * \f]
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxVarCache = elemFluxVarsCache[scvf];

        // Pressure gradient ∇_X p at the SCVF IP from Box shape gradients.
        GlobalPosition gradP(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            gradP.axpy(elemVolVars[localDof.index()].pressure(),
                       fluxVarCache.gradN(localDof.index()));

        // Deformation gradient F at SCVF IP from the coupled momentum domain.
        const Tensor F = problem.couplingManager().deformationGradientAtPoint(fvGeometry, scvf.ipGlobal());
        const Scalar J = F.determinant();

        // Permeability K(J) = K0 * exp(M*(J-1)/J).
        const Scalar K = problem.spatialParams().permeability(element, J);
        const Scalar mu_f = problem.spatialParams().fluidViscosity();

        // D_darcy = J*K/μ_f * F^{-1} * F^{-T} (pulled-back Darcy tensor)
        Tensor Finv = F; Finv.invert();
        const Tensor FinvT = transpose(Finv);
        Tensor Ddarcy = Finv;
        Ddarcy.rightmultiply(FinvT);
        Ddarcy *= J * K / mu_f;

        NumEqVector flux(0.0);
        flux[GripperMassIndices::massEqIdx] =
            -vtmv(scvf.unitOuterNormal(), Ddarcy, gradP) * scvf.area();
        return flux;
    }
};

} // end namespace Dumux

// ============================================================
// 3. Property registrations
// ============================================================

namespace Dumux::Properties::TTag {
//! \brief Type tag for the gripper momentum model (displacement, PQ1Bubble).
struct GripperMomentumModel  { using InheritsFrom = std::tuple<ModelProperties>; };
//! \brief Type tag for the gripper mass model (fluid pressure + solid bulk pressure, Box/P1).
struct GripperMassModel { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace

namespace Dumux::Properties {

// --- Momentum model traits ---
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::GripperMomentumModel>
{
    using type = GripperMomentumModelTraits<
        GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::GripperMomentumModel>
{ using type = GripperMomentumLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::GripperMomentumModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = GripperMomentumVolumeVariablesTraits<PV, MT>;
    using type = GripperMomentumVolumeVariables<Traits>;
};

// --- Pressure model traits ---
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::GripperMassModel>
{ using type = GripperMassModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::GripperMassModel>
{ using type = GripperMassLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::GripperMassModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = GripperMassVolumeVariablesTraits<PV, MT>;
    using type = GripperMassVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif // DUMUX_GRIPPER_MODEL_HH
