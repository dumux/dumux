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
#include <dumux/material/fluidmatrixinteractions/frictionlaws/swefrictionlaws.hh>

namespace Dumux
{

/*!
 * \ingroup SweModel
 * \brief Volume variables for the shallow water equations model.
 */
template <class TypeTag>
class SweVolumeVariables
{

    //using ParentType = SweVolumeVariables<Traits>;
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

   enum {
        // indices for primary variables
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        waterdepthIdx = Indices::waterdepthIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

public:
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv &scv)
    {
        priVars_ = extractDofPriVars(elemSol, scv);
        auto h = priVars_[waterdepthIdx];
        auto u = priVars_[velocityXIdx];
        auto v = priVars_[velocityYIdx];
        auto ks =  problem.spatialParams().ks(element, scv, elemSol);
        grav_ =  problem.spatialParams().grav();
        auto frictionlaw =  problem.spatialParams().frictionlaw();

        //
        auto ustar_h = computeUstarH(ks,h,grav_,frictionlaw);
        frictionH_ = 0; // calculate from primary variables and spatial params
        frictionUstarH_ = 0; // calculate from primary variables and spatial params
        bottom_ = problem.spatialParams().bottom(element,scv);
        auto ksH_ = computeKsH(ks,frictionlaw);

    }


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
     * \brief Return the friction u_starh inside the sub-control volume.
     *
     */
    Scalar frictionUstarH() const
    {
        return frictionUstarH_;
    }

     /*!
     * \brief Return the friction u_starh inside the sub-control volume.
     *
     */
    Scalar getBottom() const
    {
        return bottom_;
    }

    /*!
     * \brief Return the friction height.
     *
     */
    Scalar getKsH() const
    {
      return ksH_;
    }

    /*!
     * \brief Return the gravity constant.
     *
     */
    Scalar getGravity() const
    {
      return grav_;
    }

private:
    PrimaryVariables priVars_;
    Scalar frictionH_;
    Scalar bottom_;
    Scalar frictionUstarH_;
    Scalar ksH_;
    Scalar grav_;
};

}

#endif
