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
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the matrix and fracture problem.
 */

#ifndef DUMUX_FRACTURE_TEST_SPATIAL_PARAMS_HH
#define DUMUX_FRACTURE_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the matrix and fracture problem.
 */
template<class GridGeometry, class Scalar>
class MatrixFractureSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar, MatrixFractureSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = MatrixFractureSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    MatrixFractureSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                                const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity", 1.0);
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
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
