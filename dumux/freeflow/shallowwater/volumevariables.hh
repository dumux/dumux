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
 * \ingroup ShallowWaterModel
 * \copydoc Dumux::ShallowWaterVolumeVariables
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_VOLUME_VARIABLES_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_VOLUME_VARIABLES_HH

namespace Dumux {

/*!
 * \ingroup Shallow water model
 * \brief Volume variables for the shallow water equations model.
 */
template <class Traits>
class ShallowWaterVolumeVariables
{
    using Indices =  typename Traits::ModelTraits::Indices;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;

    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv &scv)
    {

        priVars_ = elemSol[scv.localDofIndex()];
        bedSurface_ = problem.spatialParams().bedSurface(element,scv);
    }

     /*!
     * \brief Return the extrusion factor (dummy variable).
     *
     */
    Scalar extrusionFactor() const
    { return 1.0; }

     /*!
     * \brief Return water detph h inside the sub-control volume.
     *
     */
    Scalar waterDepth() const
    {
        return priVars_[Indices::waterdepthIdx];
    }

    /*!
     * \brief Return water velocity component inside the sub-control volume.
     *
     * \param directionIndex index of the direction staring at x = 0
     */
    Scalar velocity(int directionIndex) const
    {

        return priVars_[Indices::velocityOffset + directionIndex];
    }

    /*!
     * \brief Return the bed surface inside the sub-control volume.
     *
     */
    Scalar bedSurface() const
    {
        return bedSurface_;
    }

private:
    PrimaryVariables priVars_;
    Scalar bedSurface_;
};

} // end namespace Dumux

#endif
