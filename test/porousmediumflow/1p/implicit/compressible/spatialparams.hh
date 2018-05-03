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
 * \brief The spatial params the incompressible test
 */
#ifndef DUMUX_COMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_COMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        compressible 1p model
 */
template<class TypeTag>
class OnePTestSpatialParams
: public FVSpatialParamsOneP<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                             typename GET_PROP_TYPE(TypeTag, Scalar),
                             OnePTestSpatialParams<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, OnePTestSpatialParams<TypeTag>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;
    OnePTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");

        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensUpperRight");
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
        if (isInLens_(scv.dofPosition()))
            return permeabilityLens_;
        else
            return permeability_;
    }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element &element,
                        const SubControlVolume &scv,
                        const ElementSolution &elemSol) const
    { return 0.4; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
