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
 * \ingroup KOmegaModel
 * \copydoc Dumux::KOmegaIndices
 */
#ifndef DUMUX_KOMEGA_MODELS_HH
#define DUMUX_KOMEGA_MODELS_HH

namespace Dumux {

/*!
 * \ingroup KOmegaModel
 * \brief The available eddy viscosity models
 *
 * \todo update
 * The following models are available:
 *  -# original model developed in \cite Wilcox08
 *  -# original model developed in \cite Wilcox88
 *
 * A good overview and additional models are given in \cite Patel1985a and online at NASA Turbulence Reference (https://turbmodels.larc.nasa.gov/)
 */
class KOmegaModels
{
public:
    static constexpr int wilcox08 = 0;
    static constexpr int wilcox88 = 1;
};

} // end namespace Dumux

#endif
