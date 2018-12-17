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
 * \ingroup OneTests
 * \brief Definition of the spatial parameters for the tissue problem
 */
#ifndef DUMUX_TISSUE_SPATIAL_PARAMS_HH
#define DUMUX_TISSUE_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Definition of the spatial parameters for the tissue problem
 */
template<class FVGridGeometry, class Scalar>
class TissueSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar, TissueSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = TissueSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
public:
    // export permeability type
    using PermeabilityType = Scalar;

    TissueSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        permeability_ = getParam<Scalar>("SpatialParams.PermeabilityTissue");
        porosity_ = 1.0;
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol the element solution vector
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The current finite element
     * \param scv The sub control volume
     * \param elemSol The current element solution vector
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

private:
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
