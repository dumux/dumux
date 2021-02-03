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
 * \ingroup Assembly
 * \brief Helper alias to select a local assembler based on the discretization scheme.
 */
#ifndef DUMUX_MD_LOCAL_ASSEMBLER_HH
#define DUMUX_MD_LOCAL_ASSEMBLER_HH

#include <dumux/discretization/method.hh>
#include "boxlocalassembler.hh"

namespace Dumux {
namespace Impl {

    template<std::size_t id, class Assembler, DiscretizationMethod dm, DiffMethod diffMethod>
    struct MDLocalAssemblerChooser;

    template<std::size_t id, class Assembler, DiffMethod diffMethod>
    struct MDLocalAssemblerChooser<id, Assembler, DiscretizationMethod::box, diffMethod>
    { using type = SubDomainBoxLocalAssembler<id, Assembler, diffMethod>; };

    template<std::size_t id, class Assembler, DiffMethod diffMethod>
    using MDLocalAssemblerType = typename MDLocalAssemblerChooser<id,
                                                                  Assembler,
                                                                  Assembler::template SubDomainGridGeometry<id>::discMethod,
                                                                  diffMethod>::type;

} // end namespace Immpl

/*!
 * \ingroup Assembly
 * \brief Helper alias to select the local assembler type from an assembler.
 */
template<std::size_t id, class Assembler, DiffMethod diffMethod>
using MultiDomainLocalAssembler = Impl::template MDLocalAssemblerType<id, Assembler, diffMethod>;

} // namespace Dumux

#endif
