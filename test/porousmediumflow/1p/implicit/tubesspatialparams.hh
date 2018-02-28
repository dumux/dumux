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
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
#ifndef DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH
#define DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
template<class TypeTag>
class TubesTestSpatialParams
: public FVSpatialParamsOneP<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                             typename GET_PROP_TYPE(TypeTag, Scalar),
                             TubesTestSpatialParams<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, TubesTestSpatialParams<TypeTag>>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    TubesTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
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
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const Scalar radius = this->radius(scv);
        const Scalar gamma = 2; // quadratic velocity profile (Poiseuille flow)
        return radius*radius/(2*(2+gamma));
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return 1.0; }

private:
    Scalar radius_, radiusMain_;
    static constexpr Scalar eps_ = 1e-8;
};

} // end namespace Dumux

#endif
