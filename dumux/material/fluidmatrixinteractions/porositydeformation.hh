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
 * \ingroup Fluidmatrixinteractions
 * \brief A relationship for the porosity of a porous medium under mechanical deformation.
 */
#ifndef DUMUX_MATERIAL_POROSITY_DEFORMATION_HH
#define DUMUX_MATERIAL_POROSITY_DEFORMATION_HH

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
     * \brief Calculates the porosity at a position inside an element
     * \note This assumes the primary variables to be organized such that
     *       the displacements in the different grid directions are stored
     *       in the first entries of the primary variable vector.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param element The finite element
     * \param elemSol The element solution
     * \param globalPos The global position (in the element)
     * \param refPoro The solid matrix porosity without deformation
     * \param minPoro A minimum porosity value
     * \param maxPoro A maximum porosity value
     *
     * \note \cite han2003 ( https://doi.org/10.1016/S0920-4105(03)00047-0 )
     *       provide a derivation for \f$\text{d} \phi = -(1 - \phi ) \text{d} \epsilon_v \f$.
     *       Here, \f$\epsilon_v\f$ is equal to \f$\text{div} \mathbf{u}\f$.
     *       By using an initial porosity \f$\phi_0\f$ and assuming  \f$ \epsilon_{v, 0} = 0 \f$,
     *       one obtains \f$\phi = \frac{\phi_0 - \text{div} \mathbf{u}}{1 - \text{div} \mathbf{u}}\f$,
     *       which is the formulation for the rock mechanics sign convention. Here we are
     *       using the continuum mechanics sign convention, thus, the final formula reads:
     *       \f$\phi = \frac{\phi_0 + \text{div} \mathbf{u}}{1 + \text{div} \mathbf{u}}\f$.
     */
    template< class FVGridGeom, class ElemSol >
    static Scalar evaluatePorosity(const FVGridGeom& gridGeometry,
                                   const typename FVGridGeom::GridView::template Codim<0>::Entity& element,
                                   const typename FVGridGeom::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate& globalPos,
                                   const ElemSol& elemSol,
                                   Scalar refPoro,
                                   Scalar minPoro = 0.0,
                                   Scalar maxPoro = 1.0)
    {
        // compute divergence of displacement at the given position
        Scalar divU = 0.0;
        const auto gradU = evalGradients(element, element.geometry(), gridGeometry, elemSol, globalPos);
        for (int dir = 0; dir < FVGridGeom::GridView::dimension; ++dir)
            divU += gradU[dir][dir];

        using std::clamp;
        return clamp((refPoro+divU)/(1.0+divU), minPoro, maxPoro);
    }

    /*!
     * \brief Calculates the porosity at a position inside an element
     * \note This assumes the primary variables to be organized such that
     *       the displacements in the different grid directions are stored
     *       in the first entries of the primary variable vector.
     *
     *
     * \param gridGeometry The finite volume grid geometry
     * \param element The finite element
     * \param elemSol The element solution
     * \param scv The sub-control volume
     * \param refPoro The solid matrix porosity without deformation
     * \param minPoro A minimum porosity value
     */
    template< class FVGridGeom, class ElemSol >
    static Scalar evaluatePorosity(const FVGridGeom& gridGeometry,
                                   const typename FVGridGeom::GridView::template Codim<0>::Entity& element,
                                   const typename FVGridGeom::SubControlVolume& scv,
                                   const ElemSol& elemSol,
                                   Scalar refPoro,
                                   Scalar minPoro = 0.0)
    {
        // evaluate the porosity at the scv center
        return evaluatePorosity(gridGeometry, element, scv.center(), elemSol, refPoro, minPoro);
    }
};

} // namespace Dumux

#endif
