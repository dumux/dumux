// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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

namespace Dumux {

/*!
 * \brief The available discretization methods in Dumux
 * \ingroup Discretization
 * \note Use none if specifying a discretization method is required but
 *       the class in question is not specific to a a discretization method
 *       or the classification is non-applicable
 */
enum class DiscretizationMethod
{
    none, box, cctpfa, ccmpfa, staggered, fem, fcstaggered
};


namespace DiscretizationMethods {

struct CCTpfa : public Utility::Tag<CCTpfa> {
    static std::string name() { return "cctpfa"; }
};


struct CCMpfa : public Utility::Tag<CCMpfa> {
    static std::string name() { return "ccmpfa"; }
};


struct Box : public Utility::Tag<Box> {
    static std::string name() { return "box"; }
};



struct Staggered : public Utility::Tag<Staggered> {
    static std::string name() { return "staggered"; }
};



struct FEM : public Utility::Tag<FEM> {
    static std::string name() { return "fem"; }
};



struct FCStaggered : public Utility::Tag<FCStaggered> {
    static std::string name() { return "fcstaggered"; }
};



struct None : public Utility::Tag<None> {
    static std::string name() { return "none"; }
};


inline constexpr CCTpfa cctpfa{};
inline constexpr CCMpfa ccmpfa{};
inline constexpr Box box{};
inline constexpr Staggered staggered{};
inline constexpr FEM fem{};
inline constexpr FCStaggered fcstaggered{};
inline constexpr None none{};

} // end namespace DiscretizationMethods
} // end namespace Dumux

#endif
