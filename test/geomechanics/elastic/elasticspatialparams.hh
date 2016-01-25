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
SET_TYPE_PROP(ElSpatialParams, SpatialParams, Dumux::ElSpatialParams<TypeTag>);

}

template<class TypeTag>
class ElSpatialParams
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dimWorld=GridView::dimensionworld,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:

    ElSpatialParams(const GridView &gridView)
    {
        // rock density
        rockDensity_ = 2650.0;
        // Lame Parameters
        lambda_ = 3e9;
        mu_ = 3e9;
        // Young's modulus
        E_ = 1e7;
        // Poisson's ration
        nu_ = 0.3;
     }

    ~ElSpatialParams()
    {}

    /*!
     * \brief Define the rock density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param element The finite element
     * \param scvIdx The local index of the sub-control volume where
     */
    const Scalar rockDensity(const Element &element,
                                        int scvIdx) const
    {
        return rockDensity_;
    }

    /*!
     * \brief Define the rock density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param globalPos The position for which the rock density should be returned
     */
    const Scalar rockDensity(const GlobalPosition &globalPos) const
    {
        return rockDensity_;
    }

    /*!
     * \brief Define the Lame parameters \f$\mathrm{[Pa]}\f$.
     *
     * \param element The finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume where
     */
    const Dune::FieldVector<Scalar,2> lameParams(const Element &element,
                                           const FVElementGeometry &fvGeometry,
                                           int scvIdx) const
    {
        // Lame parameters
        Dune::FieldVector<Scalar, 2> param;

        param[0] = lambda_;
        param[1] = mu_;

        return param;
    }

    /*!
     * \brief Define Young's modulus E \f$\mathrm{[Pa]}\f$.
     *
     * \param element The finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume where
     */
    const Scalar E(const Element &element,
                   const FVElementGeometry &fvGeometry,
                   int scvIdx) const
    {
        return E_;
    }

    /*!
     * \brief Define Poisson's ratio \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume where
     */
    const Scalar nu(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        return nu_;
    }

private:
    Scalar rockDensity_;
    Scalar lambda_;
    Scalar mu_;
    Scalar E_;
    Scalar nu_;
    static constexpr Scalar eps_ = 3e-6;
};
}
#endif
