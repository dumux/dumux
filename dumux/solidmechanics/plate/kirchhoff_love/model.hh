// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup KirchhoffLovePlate
 * \brief Kirchhoff-Love plate model
 *
 * In the Kirchhoff-Love model, the plate is very thin and the rotation angles
 * \f$ \boldsymbol{\theta} \f$ are identified with the gradient of the vertical deformation \f$ w \f$:
 * \f[ \boldsymbol{\theta} = \nabla w. \f]
 * The equilibrium equation reads
 * \f[ \nabla\cdot(\nabla\cdot\mathbf{M}) = F, \f]
 * where \f$ F \f$ is the out-of-plane load and the moment resultant for an isotropic material is
 * \f[
 *   \mathbf{M}(w) = -D\left\{(1-\nu)\nabla\nabla w + \nu\operatorname{tr}(\nabla\nabla w)\,\mathbf{I}\right\},
 * \f]
 * with the bending modulus \f$ D = Et^3/(12(1-\nu^2)) \f$, Young's modulus \f$ E \f$,
 * Poisson ratio \f$ \nu \f$, and plate thickness \f$ t \f$.
 * This is a fourth-order PDE in \f$ w \f$, which is not directly amenable to
 * standard lowest-order finite volume discretization.
 *
 * \par Mixed form using Helmholtz decomposition
 * To reduce the problem to a system of second-order equations, we define the
 * shear resultant vector \f$ \mathbf{q} := \nabla\cdot\mathbf{M}(\boldsymbol{\theta}) \f$
 * and apply the Helmholtz decomposition
 * \f[ \mathbf{q} = \nabla\varphi + \mathbf{J}\nabla\psi, \f]
 * where \f$ \varphi \f$ is the gradient (irrotational) potential,
 * \f$ \psi \f$ is the curl (solenoidal) potential, and
 * \f$ \mathbf{J} = \begin{bmatrix}0&1\\-1&0\end{bmatrix} \f$
 * rotates clockwise by 90°, so that
 * \f$ \mathbf{J}\nabla\psi = (\partial_y\psi,\,-\partial_x\psi)^T \f$.
 * Substituting into the equilibrium equation and the constraint
 * \f$ \nabla w - \boldsymbol{\theta} = \mathbf{0} \f$
 * (taking its divergence and curl respectively), the system reads
 * \cite Destuynder1988
 * \f{align}{
 *   \nabla\cdot\nabla\varphi &= F,\\
 *   -\nabla\cdot(\nabla w - \boldsymbol{\theta}) &= 0,\\
 *   \nabla\cdot(\mathbf{J}\boldsymbol{\theta}) &= 0,\\
 *   -\nabla\cdot(\mathbf{M}(\boldsymbol{\theta}) - \mathbf{I}\varphi - \mathbf{J}\psi) &= \mathbf{0}.
 * \f}
 * Equations (1)-(3) are the deformation-and-potentials sub-problem in the implemented order
 * \f$ (\varphi, w, \psi) \f$.
 * In particular, equations (1) and (2) are scalar second-order equations in
 * \f$ \varphi \f$ and \f$ w \f$, while equation (3) is a scalar constraint equation.
 * Equation (4) is a vector second-order equation for the rotation field
 * \f$ \boldsymbol{\theta} \f$.
 *
 * \par Boundary conditions
 * The tensor \f$ \mathbf{T} = \mathbf{M}(\boldsymbol{\theta}) - \varphi\mathbf{I} - \psi\mathbf{J} \f$
 * of equation (4) is divergence-free, and its boundary traction has the components
 * \f[
 *   \mathbf{n}\cdot\mathbf{T}\mathbf{n} = M_{nn} - \varphi,\qquad
 *   \mathbf{s}\cdot\mathbf{T}\mathbf{n} = M_{ns} - \psi,\qquad
 *   \mathbf{s} = \mathbf{J}\mathbf{n},
 * \f]
 * where \f$ \mathbf{n} \f$ is the outward unit normal. The decomposition gives
 * \f$ \partial_n\varphi = q_n + \partial_s\psi \f$. The outward fluxes of equations
 * (1)-(3), in their implemented order, are
 * \f[
 *   \left(\partial_n\varphi,\;
 *   -\mathbf{n}\cdot(\nabla w-\boldsymbol{\theta}),\;
 *   -\mathbf{s}\cdot\boldsymbol{\theta}\right).
 * \f]
 * The outward flux of equation (4) is \f$ -\mathbf{T}\mathbf{n} \f$.
 *
 * A **clamped** edge prescribes \f$ w = 0 \f$ and \f$ \boldsymbol{\theta} = \mathbf{0} \f$.
 * Prescribing \f$ \varphi = 0 \f$ on one clamped boundary component fixes its gauge.
 * The boundary values of \f$ w \f$ and \f$ \varphi \f$ replace equations (1) and (2),
 * so the transverse flux is a support reaction. If there are several entirely clamped
 * boundary components, \f$ \varphi \f$ on each additional component must be one unknown
 * constant. Its equation is the vanishing sum of compatibility residuals over that
 * component's nodes. This retains the net compatibility condition lost when equation
 * (2) is replaced at those nodes.
 *
 * A straight **simply supported** edge prescribes \f$ w = 0 \f$ in place of equation (1)
 * and \f$ \boldsymbol{\theta}\cdot\mathbf{s} = \partial_s w = 0 \f$.
 * The normal traction \f$ \mathbf{n}\cdot\mathbf{T}\mathbf{n} = -\varphi \f$
 * imposes \f$ M_{nn} = 0 \f$, while the tangential traction is a reaction.
 * Prescribing \f$ \varphi = 0 \f$ in place of equation (2) fixes the gauge on the edge.
 * Leaving the tangential rotation free and imposing the full free-edge traction instead
 * would also impose \f$ \psi = M_{ns} \f$; at supported corners this can suppress the
 * twisting-moment jump that supplies the corner reaction.
 *
 * A **free** edge prescribes the traction of the rotation sub-problem as
 * \f$ \mathbf{T}\mathbf{n} = -\varphi\,\mathbf{n} \f$, which imposes \f$ M_{nn} = 0 \f$
 * and identifies \f$ \psi \f$ with the twisting moment \f$ M_{ns} \f$.
 * The rotation Neumann value is therefore \f$ +\varphi\mathbf{n} \f$.
 * The flux of equation (1) is the Kirchhoff effective shear,
 * \f$ \partial_n\varphi = q_n + \partial_s M_{ns} = V_n \f$, which vanishes on an unloaded
 * edge. Equation (2) has zero boundary flux, while equation (3) retains its solution-dependent
 * flux \f$ -\mathbf{s}\cdot\boldsymbol{\theta} \f$; setting this flux to zero would constrain
 * the tangential rotation. At a convex corner where two free edges meet, continuity of
 * \f$ \psi \f$ enforces matching twisting moments, the corner condition without a point load.
 *
 * If the tangential rotation is prescribed on the entire boundary, a constant shift of
 * \f$ \psi \f$ leaves the equations unchanged, so one interior value must be fixed.
 * A free edge fixes this constant through \f$ \psi = M_{ns} \f$, and no additional
 * constraint on \f$ \psi \f$ is imposed. For an entirely free boundary, \f$ \varphi \f$
 * instead needs one fixed value. This fixes the potential gauge; the affine rigid-body
 * deflection modes remain and require compatible loads and separate constraints.
 *
 * \par Primary variables
 * The deformation sub-problem has three primary variables per DOF:
 * - shear gradient potential \f$ \varphi \f$
 * - vertical deformation \f$ w \f$
 * - shear curl potential \f$ \psi \f$
 *
 * The model describes static equilibrium.
 *
 * The rotation sub-problem has two primary variables per DOF:
 * - rotation component \f$ \theta_x \f$
 * - rotation component \f$ \theta_y \f$
 */

#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_MODEL_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

template<class PV, class MT>
struct KirchhoffLovePlateVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

struct KirchhoffLovePlateIndices
{
    static constexpr int shearGradPotentialIdx = 0;
    static constexpr int verticalDeformationIdx = 1;
    static constexpr int shearCurlPotentialIdx = 2;

    static constexpr int shearGradPotentialEqIdx = 0;
    static constexpr int deformationEqIdx = 1;
    static constexpr int shearCurlPotentialEqIdx = 2;
};

/*!
 * \ingroup KirchhoffLovePlate
 * \brief KirchhoffLovePlateTraits
 */
struct KirchhoffLovePlateTraits
{
    using Indices = KirchhoffLovePlateIndices;
    static constexpr int numEq() { return 3; }
};

struct KirchhoffLovePlateRotationIndices
{
    static constexpr int rotation0Idx = 0;
    static constexpr int rotation1Idx = 1;

    static constexpr int momentEq0Idx = 0;
    static constexpr int momentEq1Idx = 1;
};

/*!
 * \ingroup KirchhoffLovePlate
 * \brief KirchhoffLovePlateRotationModelTraits
 */
struct KirchhoffLovePlateRotationModelTraits
{
    using Indices = KirchhoffLovePlateRotationIndices;
    static constexpr int numEq() { return 2; }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {
struct KirchhoffLovePlateDeformation { using InheritsFrom = std::tuple<ModelProperties>; };
struct KirchhoffLovePlateRotation { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::KirchhoffLovePlateDeformation>
{ using type = KirchhoffLovePlateTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::KirchhoffLovePlateDeformation>
{ using type = KirchhoffLovePlateLocalResidualDeformation<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::KirchhoffLovePlateDeformation>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = KirchhoffLovePlateVolumeVariablesTraits<PV, MT>;
    using type = KirchhoffLovePlateDeformationVolumeVariables<Traits>;
};

////
// Rotation model
////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::KirchhoffLovePlateRotation>
{ using type = KirchhoffLovePlateRotationModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::KirchhoffLovePlateRotation>
{ using type = KirchhoffLovePlateLocalResidualRotation<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::KirchhoffLovePlateRotation>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = KirchhoffLovePlateVolumeVariablesTraits<PV, MT>;
    using type = KirchhoffLovePlateRotationVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
