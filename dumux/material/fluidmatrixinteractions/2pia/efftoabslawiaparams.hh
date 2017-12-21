// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Fluidmatrixinteractions
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws -- in this case the interfacial area surfaces --
 *        from effective to absolute saturations.
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_IA_PARAMS_HH
#define DUMUX_EFF_TO_ABS_LAW_IA_PARAMS_HH

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws -- in this case the interfacial area surfaces --
 *        from effective to absolute saturations.
 */
template <class EffIALawParamsT>
class EffToAbsLawIAParams : public EffIALawParamsT
{
    using EffIALawParams = EffIALawParamsT;
public:
    using Scalar = typename EffIALawParams::Scalar;

    EffToAbsLawIAParams()
        : EffIALawParams()
    {}
};

}

#endif
