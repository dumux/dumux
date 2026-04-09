// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MindlinReissnerPlate
 * \brief Mindlin-Reissner plate model
 *
 * Mindlin-Reissner plate model for small strains and small rotations.
 * The stress resultants are defined by integration over the plate thickness \f$ [-\tfrac{t}{2}, \tfrac{t}{2}] \f$:
 * \f[
 *   N_{\alpha\beta} := \int_{-t/2}^{t/2}\sigma_{\alpha\beta}\,dx_3,\quad
 *   M_{\alpha\beta} := \int_{-t/2}^{t/2} x_3\,\sigma_{\alpha\beta}\,dx_3,\quad
 *   Q_\alpha := \int_{-t/2}^{t/2}\sigma_{\alpha 3}\,dx_3.
 * \f]
 * Ignoring in-plane stresses, the equilibrium equations are
 * \f{align}{
 *   \nabla\cdot\mathbf{Q} &= F,\\
 *   -\nabla\cdot\mathbf{M} + \mathbf{Q} &= \mathbf{0},
 * \f}
 * where \f$ F \f$ is the out-of-plane load. Using a plane-stress constitutive law,
 * the moment resultants for an isotropic material are
 * \f[
 *   \mathbf{M}(\boldsymbol{\theta}) = -D\left\{(1-\nu)\,\boldsymbol{\varepsilon}(\boldsymbol{\theta})
 *     + \nu\operatorname{tr}(\boldsymbol{\varepsilon}(\boldsymbol{\theta}))\,\mathbf{I}\right\},
 * \f]
 * with \f$ \boldsymbol{\varepsilon}(\boldsymbol{\theta}) = \frac{1}{2}(\nabla\boldsymbol{\theta}+(\nabla\boldsymbol{\theta})^T) \f$,
 * bending modulus \f$ D = Et^3/(12(1-\nu^2)) \f$, and the shear resultants
 * \f[
 *   \mathbf{Q}(\boldsymbol{\theta}, w) = \kappa G t\,(\nabla w - \boldsymbol{\theta}),
 * \f]
 * where \f$ G = E/(2(1+\nu)) \f$ is the shear modulus and \f$ \kappa \f$ the shear correction factor.
 *
 * \par Mixed form using potentials
 * The shear resultants admit the Helmholtz decomposition
 * \f$ \mathbf{Q} = \nabla\varphi + \mathbf{J}\nabla\psi \f$
 * in terms of scalar potentials \f$ \varphi \f$ and \f$ \psi \f$.
 * With the rotation matrix \f$ \mathbf{J} = \begin{bmatrix}0&1\\-1&0\end{bmatrix} \f$, the system becomes
 * \f{align}{
 *   \nabla\cdot\nabla\varphi &= F,\\
 *   -\nabla\cdot(\nabla w - \boldsymbol{\theta}) &= -\nabla\cdot((\kappa Gt)^{-1}\nabla\varphi),\\
 *   -\nabla\cdot(\mathbf{J}\boldsymbol{\theta}) &= \nabla\cdot((\kappa Gt)^{-1}\mathbf{J} \mathbf{J}\nabla\psi),\\
 *   -\nabla\cdot(\mathbf{M}(\boldsymbol{\theta}) - \mathbf{I}\varphi - \mathbf{J}\psi) &= \mathbf{0}.
 * \f}
 * Equations (1)-(3) are the deformation-and-potentials sub-problem according to the implemented order
 * for the primary variables \f$ (\varphi, w, \psi) \f$. All three are scalar second-order equations
 * in \f$ \varphi \f$, \f$ w \f$, and \f$ \psi \f$, respectively.
 * Equation (4) is a vector second-order equation for the rotation field \f$ \boldsymbol{\theta} \f$.
 *
 * \par Primary variables
 * The deformation sub-problem has three primary variables per DOF:
 * - shear gradient potential \f$ \varphi \f$
 * - vertical deformation \f$ w \f$
 * - shear curl potential \f$ \psi \f$
 *
 * The rotation sub-problem has two primary variables per DOF:
 * - rotation component \f$ \theta_x \f$
 * - rotation component \f$ \theta_y \f$
 */

#ifndef DUMUX_MINDLIN_REISSNER_PLATE_MODEL_HH
#define DUMUX_MINDLIN_REISSNER_PLATE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/solidmechanics/plate/kirchhoff_love/model.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/volumevariables.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/couplingmanager.hh>

#include "localresidual.hh"

namespace Dumux {

template<class MDTraits>
class MindlinReissnerPlateCouplingManager
: public KirchhoffLovePlateCouplingManager<MDTraits>
{
    using ParentType = KirchhoffLovePlateCouplingManager<MDTraits>;
public:
    using ParentType::ParentType;
};

using MindlinReissnerPlateTraits = KirchhoffLovePlateTraits;
using MindlinReissnerPlateRotationModelTraits = KirchhoffLovePlateRotationModelTraits;

template<class PV, class MT>
struct MindlinReissnerPlateVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

template<class T>
using MindlinReissnerPlateDeformationVolumeVariables
    = KirchhoffLovePlateDeformationVolumeVariables<T>;

template<class T>
using MindlinReissnerPlateRotationVolumeVariables
    = KirchhoffLovePlateRotationVolumeVariables<T>;

} // end namespace Dumux

namespace Dumux::Properties::TTag {
struct MindlinReissnerPlateDeformation { using InheritsFrom = std::tuple<ModelProperties>; };
struct MindlinReissnerPlateRotation { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::MindlinReissnerPlateDeformation>
{ using type = MindlinReissnerPlateTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::MindlinReissnerPlateDeformation>
{ using type = MindlinReissnerPlateLocalResidualDeformation<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::MindlinReissnerPlateDeformation>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = MindlinReissnerPlateVolumeVariablesTraits<PV, MT>;
    using type = MindlinReissnerPlateDeformationVolumeVariables<Traits>;
};

////
// Rotation model
////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::MindlinReissnerPlateRotation>
{ using type = MindlinReissnerPlateRotationModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::MindlinReissnerPlateRotation>
{ using type = MindlinReissnerPlateLocalResidualRotation<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::MindlinReissnerPlateRotation>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = MindlinReissnerPlateVolumeVariablesTraits<PV, MT>;
    using type = MindlinReissnerPlateRotationVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
