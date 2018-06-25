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
 * \ingroup Fluidmatrixinteractions
 * \brief A relationship for the porosity of a porous medium under mechanical deformation.
 */
#ifndef DUMUX_POROSITY_DEFORMATION_HH
#define DUMUX_POROSITY_DEFORMATION_HH

#include <dumux/discretization/evalgradients.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A relationship for the porosity of a porous medium under mechanical deformation.
 *
 * \tparam Scalar The type used for scalar values
 */
template<class Scalar>
class PorosityDeformation
{
public:
    /*!
     * \brief Calculates the porosity in a sub-control volume
     * \note This assumes the primary variables to be organized
     *       such that the displacements in the different grid
     *       directions are stored in the first entries of the
     *       primary variable vector.
     *
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param element element
     * \param elemSol the element solution
     * \param scv sub control volume
     * \param refPoro The solid matrix porosity without deformation
     * \param minPoro A minimum porosity value
     */
    template< class FVGridGeom, class ElemSol >
    Scalar evaluatePorosity(const FVGridGeom& fvGridGeometry,
                            const typename FVGridGeom::GridView::template Codim<0>::Entity& element,
                            const typename FVGridGeom::SubControlVolume& scv,
                            const ElemSol& elemSol,
                            Scalar refPoro,
                            Scalar minPoro = 0.0) const
    {
        // compute divergence of diplacement for this scv
        Scalar divU = 0.0;
        const auto gradU = evalGradients(element, element.geometry(), fvGridGeometry, elemSol, scv.center());
        for (int dir = 0; dir < FVGridGeom::GridView::dimension; ++dir)
            divU += gradU[dir][dir];

        using std::max;
        return max(minPoro, refPoro*(1.0+divU));
    }
};

} // namespace Dumux

#endif
