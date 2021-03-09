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
 * \brief Convenience alias to select a local assembler based on the discretization scheme.
 */
#ifndef DUMUX_LOCAL_ASSEMBLER_HH
#define DUMUX_LOCAL_ASSEMBLER_HH

#include <dumux/discretization/method.hh>
#include "fv/boxlocalassembler.hh"
#include "fv/cclocalassembler.hh"

namespace Dumux::Experimental {
namespace Detail {

    template<class LO, DiscretizationMethod dm>
    struct LocalAssemblerChooser;

    template<class LO>
    struct LocalAssemblerChooser<LO, DiscretizationMethod::box>
    { using type = BoxLocalAssembler<LO>; };

    template<class LO>
    struct LocalAssemblerChooser<LO, DiscretizationMethod::cctpfa>
    { using type = CCLocalAssembler<LO>; };

    template<class LO>
    struct LocalAssemblerChooser<LO, DiscretizationMethod::ccmpfa>
    { using type = CCLocalAssembler<LO>; };

    template<class LocalOperator>
    using LocalAssemblerType
        = typename LocalAssemblerChooser<LocalOperator,
            LocalOperator::LocalContext::ElementGridGeometry::GridGeometry::discMethod>::type;

} // end namespace Detail

/*!
 * \ingroup Assembly
 * \brief Helper alias to select the local assembler type from a local operator.
 */
template<class LocalOperator>
using LocalAssembler = Detail::LocalAssemblerType<LocalOperator>;

} // end namespace Dumux::Experimental

#endif
