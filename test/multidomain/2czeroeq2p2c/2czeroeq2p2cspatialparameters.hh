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
 * \brief Spatial parameters for the
 *        coupling of an isothermal two-component ZeroEq
 *        and an isothermal two-phase two-component Darcy model.
 */

#ifndef DUMUX_TWOCZEROEQTWOPTWOCSPATIALPARAMS_HH
#define DUMUX_TWOCZEROEQTWOPTWOCSPATIALPARAMS_HH

#include <dune/grid/io/file/vtk/common.hh>

#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoCZeroEqTwoPTwoCSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoCZeroEqTwoPTwoCSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCSpatialParams, SpatialParams, TwoCZeroEqTwoPTwoCSpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedVanGenuchten<typename GET_PROP_TYPE(TypeTag, Scalar)>>);
//               EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}


/*!
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for
 *        the coupling of an isothermal two-component ZeroEq
 *        and an isothermal two-phase two-component Darcy model.
 */
template<class TypeTag>
class TwoCZeroEqTwoPTwoCSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief Spatial parameters for the
     *        coupling of an isothermal two-component ZeroEq
     *        and an isothermal two-phase two-component Darcy model.
     *
     * \param gridView The GridView which is used by the problem
     */
    TwoCZeroEqTwoPTwoCSpatialParams(const GridView& gridView)
        : ParentType(gridView)
    {
        permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Permeability);
        porosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Porosity);

        spatialParams_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Swr));
        spatialParams_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Snr));
        spatialParams_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, VgAlpha));
        spatialParams_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, VgN));
    }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       const int scvIdx) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the parameter object for the material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
        return spatialParams_;
    }

private:
    Scalar permeability_;
    Scalar porosity_;
    MaterialLawParams spatialParams_;
};

} // end namespace Dumux

#endif // DUMUX_TWOCZEROEQTWOPTWOCSPATIALPARAMS_HH
