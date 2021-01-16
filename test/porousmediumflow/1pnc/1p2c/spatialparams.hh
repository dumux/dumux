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
 * \ingroup OnePNCTests
 * \brief Definition of the spatial parameters for the 1pnc problems.
 */

#ifndef DUMUX_1PNC_TEST_SPATIAL_PARAMS_HH
#define DUMUX_1PNC_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Definition of the spatial parameters for the 1pnc test problems.
 */
template<class GridGeometry, class Scalar>
class OnePNCTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             OnePNCTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           OnePNCTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePNCTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Defines the dispersivity.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The solution for all dofs of the element
     */
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }

private:
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
