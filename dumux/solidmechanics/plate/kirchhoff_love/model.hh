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
 * is the 90° rotation matrix (so that \f$\mathbf{J}\nabla\psi\f$ is the
 * 2D curl of \f$\psi\f$, i.e. \f$(-\partial_y\psi,\,\partial_x\psi)^T\f$).
 * Substituting into the equilibrium equation and the constraint
 * \f$ \nabla w - \boldsymbol{\theta} = \mathbf{0} \f$
 * (taking its divergence and curl respectively), the system reads
 * \cite Destuynder1988
 * \f{align}{
 *   \nabla\cdot\nabla\varphi &= F,\\
 *   -\nabla\cdot(\nabla w - \boldsymbol{\theta}) &= 0,\\
 *   -\nabla\cdot(\mathbf{J}\boldsymbol{\theta}) &= 0,\\
 *   -\nabla\cdot(\mathbf{M}(\boldsymbol{\theta}) - \mathbf{I}\varphi - \mathbf{J}\psi) &= \mathbf{0}.
 * \f}
 * Equations (1)-(3) are the deformation-and-potentials sub-problem in the implemented order
 * \f$ (\varphi, w, \psi) \f$.
 * In particular, equations (1) and (2) are scalar second-order equations in
 * \f$ \varphi \f$ and \f$ w \f$, while equation (3) is a scalar constraint equation.
 * Equation (4) is a vector second-order equation for the rotation field
 * \f$ \boldsymbol{\theta} \f$.
 *
 * \par Primary variables
 * The deformation sub-problem has three primary variables per DOF:
 * - shear gradient potential \f$ \varphi \f$
 * - vertical deformation \f$ w \f$
 * - shear curl potential \f$ \psi \f$
 *
 * The current implementation solves only the static (equilibrium) problem.
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
