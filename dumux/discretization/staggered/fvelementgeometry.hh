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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFVElementGeometry
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH

#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This locally builds up the sub control volumes and sub control volume faces
 *        for each element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableFVGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableFVGridGeometryCache>
class StaggeredFVElementGeometry : public CCTpfaFVElementGeometry<GG, enableFVGridGeometryCache>
{
    using ParentType = CCTpfaFVElementGeometry<GG, enableFVGridGeometryCache>;
    using IndexType = typename GG::GridView::IndexSet::IndexType;
public:
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    using ParentType::ParentType;

    //! Constructor getting a auxiliary cell center of face specific FvGridGeometry type.
    //! Needed for the multi-domain framework.
    template<class CellCenterOrFaceFVGridGeometry>
    StaggeredFVElementGeometry(const CellCenterOrFaceFVGridGeometry& fvGridGeometry)
    : ParentType(fvGridGeometry.actualfvGridGeometry()) {}

    //! Get a sub control volume face with an element index and a local scvf index
    const SubControlVolumeFace& scvf(IndexType eIdx, IndexType localScvfIdx) const
    {
        return this->fvGridGeometry().scvf(eIdx, localScvfIdx);
    }
};

} // end namespace

#endif
