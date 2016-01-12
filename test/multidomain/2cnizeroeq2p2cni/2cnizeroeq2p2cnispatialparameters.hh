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
 *        coupling of a non-isothermal two-component ZeroEq
 *        and a non-isothermal two-phase two-component Darcy model.
 */
#ifndef DUMUX_TWOCNIZEROEQTWOPTWOCNISPATIALPARAMS_HH
#define DUMUX_TWOCNIZEROEQTWOPTWOCNISPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class TwoCNIZeroEqTwoPTwoCNISpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoCNIZeroEqTwoPTwoCNISpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNISpatialParams, SpatialParams,
              TwoCNIZeroEqTwoPTwoCNISpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNISpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedVanGenuchten<typename GET_PROP_TYPE(TypeTag, Scalar)>>);
//               EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}


/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 * \brief Definition of the spatial parameters for
 *        the coupling of a non-isothermal two-component ZeroEq
 *        and a non-isothermal two-phase two-component Darcy model.
 */
template<class TypeTag>
class TwoCNIZeroEqTwoPTwoCNISpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

public:
    /*!
     * \brief Spatial parameters for the
     *        coupling of a non-isothermal two-component ZeroEq
     *        and a non-isothermal two-phase two-component Darcy model.
     *
     * \param gridView The GridView which is used by the problem
     */
    TwoCNIZeroEqTwoPTwoCNISpatialParams(const GridView& gridView)
        : ParentType(gridView)
    {
        permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Permeability);
        porosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Porosity);
        thermalConductivitySolid_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, ThermalConductivitySolid);
        alphaBJ_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, AlphaBJ);

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
    Scalar intrinsicPermeability(const Element &element,
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

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidHeatCapacity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return 790.0;
    }

    /*!
     * \brief Returns the density of the solid material (not the bulk density) \f$[kg/m^3]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidDensity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return 2650.0;
    }

    /*!
     * \brief Returns the thermal conductivity \f$[W/(m*K)]\f$ of the solid
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        return thermalConductivitySolid_;
    }

    /*!
     * \brief Evaluate the Beavers-Joseph coefficient at given position
     *
     * \param globalPos The global position
     *
     * \return Beavers-Joseph coefficient
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition &globalPos) const
    {
        return alphaBJ_;
    }

private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar thermalConductivitySolid_;
    Scalar alphaBJ_;

    MaterialLawParams spatialParams_;
};
} // end namespace

#endif // DUMUX_TWOCNIZEROEQTWOPTWOCNISPATIALPARAMS_HH
