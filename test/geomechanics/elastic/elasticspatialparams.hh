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
 * \brief Definition of the spatial parameters for the linear elasticity problem.
 */
#ifndef DUMUX_ELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_ELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/el2p/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ElasticBoxModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the linear elasticity
 *        problem.
 */

template<class TypeTag>
class ElSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ElSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ElSpatialParams, SpatialParams, ElSpatialParams<TypeTag>);
}

template<class TypeTag>
class ElSpatialParams
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    using CoordScalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dimWorld = GridView::dimensionworld;
    using LameParams = Dune::FieldVector<Scalar,2>;
    using GlobalPosition = Dune::FieldVector<CoordScalar,dimWorld>;

public:

    //! The constructor
    ElSpatialParams(const Problem& problem, const GridView &gridView) {}

    /*!
     * \brief Define the rock density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param element The finite element
     * \param ipSol The solution at the integration point
     */
    Scalar rockDensity(const Element &element,
                       const PrimaryVariables& ipSol) const
    { return 2650.0; }

    /*!
     * \brief Define the Lame parameters \f$\mathrm{[Pa]}\f$.
     *
     * \param element The finite element
     * \param ipSol The solution at the integration point
     */
    LameParams lameParams(const Element &element,
                          const PrimaryVariables& ipSol) const
    { return LameParams({3e9, 3e9}); }

    /*!
     * \brief Define Young's modulus E \f$\mathrm{[Pa]}\f$.
     *
     * \param element The finite element
     * \param ipSol The solution at the integration point
     */
    Scalar E(const Element &element,
             const PrimaryVariables& ipSol) const
    { return 1e7; }

    /*!
     * \brief Define Poisson's ratio \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param ipSol The solution at the integration point
     */
    Scalar nu(const Element &element,
              const PrimaryVariables& ipSol) const
    { return 0.3; }
};
}
#endif
