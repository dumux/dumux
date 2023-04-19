// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief The available discretization methods in Dumux
 */
#ifndef DUMUX_DISCRETIZATION_METHOD_HH
#define DUMUX_DISCRETIZATION_METHOD_HH

#include <ostream>
#include <string>

#include <dumux/common/tag.hh>

namespace Dumux::DiscretizationMethods {

/*
 * \brief Cell-centered finite volume scheme with two-point flux approximation
 */
struct CCTpfa : public Utility::Tag<CCTpfa> {
    static std::string name() { return "cctpfa"; }
};


/*
 * \brief Cell-centered finite volume scheme with multi-point flux approximation
 */
struct CCMpfa : public Utility::Tag<CCMpfa> {
    static std::string name() { return "ccmpfa"; }
};


/*
 * \brief Control-volume finite element methods
 * This is a group of discretization methods that share certain properties.
 * Therefore there is a single meta-tag parametrized in terms of the actual
 * discretization method in the group. Having a common tag allows to specialize
 * template agnostic of the underlying discretization type
 */
template<class DM>
struct CVFE : public Utility::Tag<CVFE<DM>> {
    static std::string name() { return DM::name(); }
};


#ifndef DOXYGEN
namespace Detail {

template<class DM>
struct IsCVFE : public std::false_type {};

template<class DM>
struct IsCVFE<CVFE<DM>> : public std::true_type {};

} // end namespace Detail
#endif

/*
 * \brief Template variable that is true when the discretization method DM is a CVFE schemes
 */
template<class DM>
inline constexpr bool isCVFE = Detail::IsCVFE<DM>::value;


/*
 * \brief Various control volume finite element discretization methods
 */
namespace CVFEMethods {

struct PQ1 {
    static std::string name() { return "box"; }
};

struct CR_RT {
    static std::string name() { return "fcdiamond"; }
};

struct PQ1Bubble {
    static std::string name() { return "pq1bubble"; }
};

} // end namespace CVFEMethods


/*
 * \brief Vertex-centered finite volume scheme
 * or control-volume finite element scheme based on a P1 (simplices) or Q1 (quads) basis
 */
using Box = CVFE<CVFEMethods::PQ1>;

/*
 * \brief Face-centered finite volume scheme
 * or control-volume finite element scheme based on
 * Crouzeix-Raviart (simplices) or Rannacher-Turek (quads) basis
 */
using FCDiamond = CVFE<CVFEMethods::CR_RT>;

/*
 * \brief Vertex- and cell-centered finite volume scheme
 * or control-volume finite element scheme based on
 * linear Lagrangian elements with bubble function
 */
using PQ1Bubble = CVFE<CVFEMethods::PQ1Bubble>;


/*
 * \brief Staggered-grid finite volume scheme (old)
 */
struct Staggered : public Utility::Tag<Staggered> {
    static std::string name() { return "staggered"; }
};


/*
 * \brief Finite element method
 */
struct FEM : public Utility::Tag<FEM> {
    static std::string name() { return "fem"; }
};


/*
 * \brief Staggered-grid finite volume scheme
 */
struct FCStaggered : public Utility::Tag<FCStaggered> {
    static std::string name() { return "fcstaggered"; }
};


/*
 * \brief Tag used for defaults not depending on the discretization
 * or in situations where a discretization tag is needed but none
 * can be provided (the implementation has to support this of course)
 */
struct None : public Utility::Tag<None> {
    static std::string name() { return "none"; }
};


inline constexpr CCTpfa cctpfa{};
inline constexpr CCMpfa ccmpfa{};
inline constexpr Box box{};
inline constexpr PQ1Bubble pq1bubble{};
inline constexpr Staggered staggered{};
inline constexpr FEM fem{};
inline constexpr FCStaggered fcstaggered{};
inline constexpr FCDiamond fcdiamond{};
inline constexpr None none{};

} // end namespace Dumux::DiscretizationMethods

#endif
