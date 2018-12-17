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
  * \ingroup CCMpfaDiscretization
  * \brief Data structure holding interaction volume-local information for
  *        a grid subb-control volume face embedded in it.
  */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_FACE_DATA_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_FACE_DATA_HH

#include <cassert>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief General implementation of a data structure holding interaction
 *        volume-local information for a grid sub-control volume face embedded in it.
 *
 * \tparam GridIndexType The type used for indices on the grid
 * \tparam LocalIndexType The type used for indices inside interaction volumes
 */
template< class GridIndexType, class LocalIndexType >
class InteractionVolumeLocalFaceData
{
    LocalIndexType ivLocalScvfIndex_;          //!< the iv-local scvf index this scvf maps to
    LocalIndexType ivLocalInsideScvIndex_;     //!< the iv-local index of the scvfs' inside scv
    LocalIndexType scvfLocalOutsideScvfIndex_; //!< the index of this scvf in the scvf-local outside faces
    GridIndexType gridScvfIndex_;              //!< the index of the corresponding global scvf
    bool isOutside_;                           //!< indicates if this face is an "outside" face in the iv-local system

public:
    //! Default constructor
    InteractionVolumeLocalFaceData() = default;

    //! Constructor
    InteractionVolumeLocalFaceData(LocalIndexType faceIndex,
                                   LocalIndexType scvIndex,
                                   GridIndexType gridScvfIndex)
    : ivLocalScvfIndex_(faceIndex)
    , ivLocalInsideScvIndex_(scvIndex)
    , gridScvfIndex_(gridScvfIndex)
    , isOutside_(false)
    {}

    //! Constructor for "outside" faces
    InteractionVolumeLocalFaceData(LocalIndexType faceIndex,
                                   LocalIndexType scvIndex,
                                   LocalIndexType indexInScvfOutsideFaces,
                                   GridIndexType gridScvfIndex)
    : ivLocalScvfIndex_(faceIndex)
    , ivLocalInsideScvIndex_(scvIndex)
    , scvfLocalOutsideScvfIndex_(indexInScvfOutsideFaces)
    , gridScvfIndex_(gridScvfIndex)
    , isOutside_(true)
    {}

    // Functions to return stored data
    LocalIndexType ivLocalScvfIndex() const { return ivLocalScvfIndex_; }
    LocalIndexType ivLocalInsideScvIndex() const { return ivLocalInsideScvIndex_; }
    LocalIndexType scvfLocalOutsideScvfIndex() const { assert(isOutside_); return scvfLocalOutsideScvfIndex_; }
    GridIndexType gridScvfIndex() const { return gridScvfIndex_; }
    bool isOutsideFace() const { return isOutside_; }
};

} // end namespace Dumux

#endif
