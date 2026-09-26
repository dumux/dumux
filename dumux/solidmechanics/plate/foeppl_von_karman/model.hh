// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FoepplVonKarmanPlate
 * \brief Föppl-von Kármán plate model
 *
 * The Föppl-von Kármán model extends the Kirchhoff-Love plate (@ref KirchhoffLovePlate)
 * to moderately large deflections, where the vertical deformation \f$ w \f$ is of the
 * order of the plate thickness \f$ t \f$. The rotations remain small, but the
 * membrane strain picks up the quadratic contribution of the slope,
 * \f[
 *   \boldsymbol{\varepsilon}(\mathbf{u}, w)
 *     = \tfrac{1}{2}\left(\nabla\mathbf{u} + (\nabla\mathbf{u})^T
 *                         + \nabla w\otimes\nabla w\right),
 * \f]
 * with the in-plane displacement \f$ \mathbf{u} = (u_1, u_2) \f$. For an isotropic
 * material in plane stress the in-plane stress resultant is
 * \f[
 *   \mathbf{N}(\mathbf{u}, w) = C\left\{(1-\nu)(\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}_g)
 *     + \nu\operatorname{tr}(\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}_g)\,\mathbf{I}\right\},\quad
 *   C = \frac{Et}{1-\nu^2},
 * \f]
 * where the eigenstrain \f$ \boldsymbol{\varepsilon}_g \f$ is the incompatible strain of
 * growth or thermal expansion, provided by the in-plane and deformation problems through
 * an optional method `eigenstrain(globalPos)` and zero otherwise,
 * and the moment resultant is as in the Kirchhoff-Love model,
 * \f[
 *   \mathbf{M}(w) = -D\left\{(1-\nu)\nabla\nabla w + \nu\operatorname{tr}(\nabla\nabla w)\,\mathbf{I}\right\},\quad
 *   D = \frac{Et^3}{12(1-\nu^2)}.
 * \f]
 * The two equilibrium equations are coupled in both directions,
 * \f{align}{
 *   \nabla\cdot(\nabla\cdot\mathbf{M}) + \nabla\cdot(\mathbf{N}\nabla w) &= F,\\
 *   -\nabla\cdot\mathbf{N} &= \mathbf{f},
 * \f}
 * where \f$ F \f$ is the out-of-plane load and \f$ \mathbf{f} \f$ an in-plane body force
 * (both per unit area). The stretching of the mid-surface stiffens the plate against
 * bending, and conversely the deflection stretches the mid-surface.
 *
 * \par Mixed form using Helmholtz decomposition
 * The bending part is split exactly as in the Kirchhoff-Love model: with the shear
 * resultant \f$ \mathbf{q} := \nabla\cdot\mathbf{M}(\boldsymbol{\theta}) \f$ and its
 * Helmholtz decomposition \f$ \mathbf{q} = \nabla\varphi + \mathbf{J}\nabla\psi \f$,
 * \f$ \mathbf{J} = \begin{bmatrix}0&1\\-1&0\end{bmatrix} \f$, the system reads
 * \f{align}{
 *   \nabla\cdot(\nabla\varphi + \mathbf{N}\nabla w) &= F,\\
 *   -\nabla\cdot(\nabla w - \boldsymbol{\theta}) &= 0,\\
 *   -\nabla\cdot(\mathbf{J}\boldsymbol{\theta}) &= 0,\\
 *   -\nabla\cdot(\mathbf{M}(\boldsymbol{\theta}) - \mathbf{I}\varphi - \mathbf{J}\psi) &= \mathbf{0},\\
 *   -\nabla\cdot\mathbf{N}(\mathbf{u}, w) &= \mathbf{f}.
 * \f}
 * Only equation (1) differs from the Kirchhoff-Love system, since
 * \f$ \nabla\cdot(\mathbf{J}\nabla\psi) = 0 \f$ leaves the decomposition of
 * \f$ \mathbf{q} \f$ untouched by the membrane term. Setting \f$ \mathbf{N} = \mathbf{0} \f$
 * recovers the Kirchhoff-Love model; setting \f$ D = 0 \f$ and \f$ \mathbf{N} = T\mathbf{I} \f$
 * recovers the membrane model (@ref MembranePlate).
 *
 * \par Boundary conditions
 * As in the Kirchhoff-Love model (@ref KirchhoffLovePlate), where the traction
 * \f$ \mathbf{T}\mathbf{n} = -\varphi\mathbf{n} \f$ of the rotation sub-problem is what
 * makes an edge free. The flux of equation (1) now carries the membrane traction as well and
 * is the transverse force the edge transmits,
 * \f$ \partial_n\varphi + \mathbf{n}\cdot\mathbf{N}\nabla w \f$. On a free edge its
 * membrane part vanishes with the in-plane traction \f$ \mathbf{N}\mathbf{n} = \mathbf{0} \f$,
 * so the free edge is again a zero Neumann condition on the whole deformation sub-problem,
 * together with \f$ \mathbf{T}\mathbf{n} = -\varphi\mathbf{n} \f$ on the rotations and a
 * traction-free in-plane sub-problem.
 *
 * \par Primary variables
 * The model is implemented as three coupled sub-problems.
 * The deformation sub-problem has three primary variables per DOF:
 * - shear gradient potential \f$ \varphi \f$
 * - vertical deformation \f$ w \f$
 * - shear curl potential \f$ \psi \f$
 *
 * The rotation sub-problem has two primary variables per DOF:
 * - rotation component \f$ \theta_1 \f$
 * - rotation component \f$ \theta_2 \f$
 *
 * The in-plane sub-problem has two primary variables per DOF:
 * - in-plane displacement component \f$ u_1 \f$
 * - in-plane displacement component \f$ u_2 \f$
 *
 * \par Dynamics
 * The transverse balance carries the inertia \f$ \rho h\,\ddot w \f$ of the plate when the
 * deformation problem provides a mass per area and an acceleration, the latter evaluated by a
 * time integration scheme the problem holds, such as the Newmark scheme. Rotary inertia is neglected, as in the Kirchhoff-Love kinematics. Without those two
 * methods the model is static, and the optional damping term of the deformation problem then
 * turns a time loop into a relaxation towards a stable equilibrium.
 */

#ifndef DUMUX_FOEPPL_VON_KARMAN_PLATE_MODEL_HH
#define DUMUX_FOEPPL_VON_KARMAN_PLATE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/solidmechanics/plate/kirchhoff_love/model.hh>

#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

template<class PV, class MT>
struct FoepplVonKarmanPlateVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

using FoepplVonKarmanPlateIndices = KirchhoffLovePlateIndices;
using FoepplVonKarmanPlateTraits = KirchhoffLovePlateTraits;
using FoepplVonKarmanPlateRotationIndices = KirchhoffLovePlateRotationIndices;
using FoepplVonKarmanPlateRotationModelTraits = KirchhoffLovePlateRotationModelTraits;

struct FoepplVonKarmanPlateInPlaneIndices
{
    static constexpr int displacement0Idx = 0;
    static constexpr int displacement1Idx = 1;

    static constexpr int forceEq0Idx = 0;
    static constexpr int forceEq1Idx = 1;
};

/*!
 * \ingroup FoepplVonKarmanPlate
 * \brief FoepplVonKarmanPlateInPlaneModelTraits
 */
struct FoepplVonKarmanPlateInPlaneModelTraits
{
    using Indices = FoepplVonKarmanPlateInPlaneIndices;
    static constexpr int numEq() { return 2; }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {
struct FoepplVonKarmanPlateDeformation { using InheritsFrom = std::tuple<ModelProperties>; };
struct FoepplVonKarmanPlateRotation { using InheritsFrom = std::tuple<KirchhoffLovePlateRotation>; };
struct FoepplVonKarmanPlateInPlane { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::FoepplVonKarmanPlateDeformation>
{ using type = FoepplVonKarmanPlateTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FoepplVonKarmanPlateDeformation>
{ using type = FoepplVonKarmanPlateLocalResidualDeformation<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::FoepplVonKarmanPlateDeformation>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = FoepplVonKarmanPlateVolumeVariablesTraits<PV, MT>;
    using type = FoepplVonKarmanPlateDeformationVolumeVariables<Traits>;
};

////
// In-plane model
////

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::FoepplVonKarmanPlateInPlane>
{ using type = FoepplVonKarmanPlateInPlaneModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FoepplVonKarmanPlateInPlane>
{ using type = FoepplVonKarmanPlateLocalResidualInPlane<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::FoepplVonKarmanPlateInPlane>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = FoepplVonKarmanPlateVolumeVariablesTraits<PV, MT>;
    using type = FoepplVonKarmanPlateInPlaneVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
