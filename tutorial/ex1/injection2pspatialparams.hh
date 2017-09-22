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
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotmateriallaw.hh>

#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class InjectionSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InjectionSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InjectionSpatialParams, SpatialParams, InjectionSpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(InjectionSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */
template<class TypeTag>
class InjectionSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InjectionSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        aquiferHeightFromBottom_ = 30.0;

        // intrinsic permeabilities
        aquitardK_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PermeabilityAquitard);
        aquiferK_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PermeabilityAquifer);

        // porosities
        aquitardPorosity_ = 0.2;
        aquiferPorosity_ = 0.4;

        // residual saturations
        aquitardMaterialParams_.setSwr(0.2);
        aquitardMaterialParams_.setSnr(0.0);
        aquiferMaterialParams_.setSwr(0.2);
        aquiferMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        aquitardMaterialParams_.setPe(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.EntryPressureAquitard));
        aquiferMaterialParams_.setPe(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.EntryPressureAquifer));
        aquitardMaterialParams_.setLambda(2.0);
        aquiferMaterialParams_.setLambda(2.0);

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
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        if (isInAquitard_(globalPos))
            return aquitardK_;
        return aquiferK_;
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
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        if (isInAquitard_(globalPos))
            return aquitardPorosity_;
        return aquiferPorosity_;
    }


    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        if (isInAquitard_(globalPos))
            return aquitardMaterialParams_;
        return aquiferMaterialParams_;
    }

    /*!
     *  These parameters are only needed for nonisothermal models. Comment them in if you want to implement the 2pni model.
     */

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
//     Scalar solidHeatCapacity(const Element &element,
//                              const FVElementGeometry &fvGeometry,
//                              const int scvIdx) const
//     {
//         return 790; // specific heat capacity of granite [J / (kg K)]
//     }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
//     Scalar solidDensity(const Element &element,
//                         const FVElementGeometry &fvGeometry,
//                         const int scvIdx) const
//     {
//         return 2700; // density of granite [kg/m^3]
//     }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the solid
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
//     Scalar solidThermalConductivity(const Element &element,
//                                     const FVElementGeometry &fvGeometry,
//                                     const int scvIdx) const
//     {
//         return lambdaSolid_;
//     }

private:

    static constexpr Scalar eps_ = 1e-6;

    bool isInAquitard_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > aquiferHeightFromBottom_ + eps_; }

    Scalar aquitardK_;
    Scalar aquiferK_;
    Scalar aquiferHeightFromBottom_;


    Scalar aquitardPorosity_;
    Scalar aquiferPorosity_;

    MaterialLawParams aquitardMaterialParams_;
    MaterialLawParams aquiferMaterialParams_;
};

}

#endif
