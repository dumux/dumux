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
 * \ingroup Adaptive
 * \brief Interface to be used by classes transferring grid data on adpative grids
 */
#ifndef DUMUX_ADAPTIVE_GRIDDATATRANSFER_HH
#define DUMUX_ADAPTIVE_GRIDDATATRANSFER_HH

namespace Dumux {

/*!
 * \ingroup Adaptive
 * \brief Interface to be used by classes transferring grid data on adpative grids
 */
class GridDataTransfer
{
public:
    //! pure virtual base class needs virtual destructor
    virtual ~GridDataTransfer() = default;

    //! store user data before grid adaption
    virtual void store() = 0;

    //! store user data after adaption
    virtual void reconstruct() = 0;
};

} // end namespace Dumux

#endif
