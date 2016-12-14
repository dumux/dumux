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
 *
  * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
#ifndef DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH
#define DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
template<class TypeTag>
class TubesTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

public:
    TubesTestSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {
        radius_ = 1.0;

        using std::sqrt;
        radiusMain_ = sqrt(sqrt(4.0/sqrt(3.0)));
    }

    /*!
     * \brief Return the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param scv The sub control volume
     */
    Scalar radius(const SubControlVolume &scv) const
    {
        if(scv.center()[2] > 0.5 - eps_)
            return radiusMain_;
        else
            return radius_;
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability.
     *
     * \param scv The sub control volume
     * \param volVars The volume variables in the sub control volumes
     * \return the intrinsic permeability
     */
    Scalar intrinsicPermeability (const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
        const Scalar radius = this->radius(scv);
        const Scalar gamma = 2; // quadratic velocity profile (Poiseuille flow)
        return radius*radius/(2*(2+gamma));
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param scv The sub control volume
     * \return the porosity
     */
    Scalar porosity(const SubControlVolume &scv) const
    { return 1.0; }

private:
    Scalar radius_, radiusMain_;
    static constexpr Scalar eps_ = 1e-8;
};

} // end namespace Dumux

#endif
