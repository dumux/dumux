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
 * \ingroup SweModel
 *
 * \copydoc Dumux::SweVolumeVariables
 */
#ifndef DUMUX_SWE_VOLUME_VARIABLES_HH
#define DUMUX_SWE_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

// forward declaration
template <class TypeTag>
class SweVolumeVariablesImplementation;

/*!
 * \ingroup SweModel
 * \brief Volume variables for the shallow water equations
 *
 */
template <class TypeTag>
using SweVolumeVariables = SweVolumeVariablesImplementation<TypeTag>;
/*!
 * \ingroup SweModel
 * \brief Volume variables for the shallow water equations model.
 */
template <class TypeTag>
class SweVolumeVariablesImplementation
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

   enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices for primary variables
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        waterdepthIdx = Indices::waterdepthIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

public:

    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        priVars_ = extractDofPriVars(elemSol, scv);
        frictionH_ = 0; // calculate from primary variables and spatial params
        frictionUstarH_ = 0; // calculate from primary variables and spatial params
    }


    /*!
     * \brief Returns the primary variables at the dof associated with a given scv.
     */
    static const auto& extractDofPriVars(const ElementSolutionVector& elemSol,
                                         const SubControlVolume& scv)
    { return elemSol[0]; }


    Scalar extrusionFactor() const
    { return 1.0; }

    /*!
     * \brief Return water detph h inside the sub-control volume.
     *
     */
    Scalar getH() const
    {
        return priVars_[waterdepthIdx];
    }

    /*!
     * \brief Return water velocity u inside the sub-control volume.
     *
     */
    Scalar getU() const
    {
        return priVars_[velocityXIdx];
    }
    /*!
     * \brief Return water velocity v inside the sub-control volume.
     *
     */
    Scalar getV() const
    {
        return priVars_[velocityYIdx];
    }

    /*!
     * \brief Return the friction height h inside the sub-control volume.
     *
     */
    Scalar frictionH() const
    {
        return frictionH_;
    }


    /*!
     * \brief Return the friction u_starh inside the sub-control volume.
     *
     */
    Scalar frictionUstarH() const
    {
        return frictionUstarH_;
    }

private:
    PrimaryVariables priVars_;
    Scalar frictionH_;
    Scalar frictionUstarH_;
};

}

#endif
