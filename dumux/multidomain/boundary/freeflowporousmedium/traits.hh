// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Traits for the free-flow/porous-medium-flow coupling.
 */

#ifndef DUMUX_MD_FREEFLOW_POROUSMEDIUM_TRAITS_HH
#define DUMUX_MD_FREEFLOW_POROUSMEDIUM_TRAITS_HH

#include <type_traits>
#include <dumux/flux/fickslaw_fwd.hh>

namespace Dumux::FreeFlowPorousMediumCoupling {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 *        Specialization for the case of using an adapter only for the free-flow model.
 * \tparam FFFS The free-flow fluidsystem
 * \tparam PMFS The porous-medium flow fluidsystem
 */
template<class FFFS, class PMFS>
struct IsSameFluidSystem
{
    static_assert(FFFS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<typename FFFS::MultiPhaseFluidSystem, PMFS>::value;
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 * \tparam FS The fluidsystem
 */
template<class FS>
struct IsSameFluidSystem<FS, FS>
{
    static_assert(FS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<FS, FS>::value; // always true
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs indicates that Fick's law is not used for diffusion.
 * \tparam DiffLaw The diffusion law.
 */
template<class DiffLaw>
struct IsFicksLaw : public std::false_type {};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs indicates that Fick's law is used for diffusion.
 * \tparam DiffLaw The diffusion law.
 */
template<class T>
struct IsFicksLaw<Dumux::FicksLaw<T>> : public std::true_type {};

} // end namespace Dumux::FreeFlowPorousMediumCoupling

#endif
