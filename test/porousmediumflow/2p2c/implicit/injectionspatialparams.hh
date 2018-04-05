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
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

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
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */
template<class TypeTag>
class InjectionSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static constexpr int dimWorld = GridView::dimensionworld;

    using CoordScalar = typename Grid::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    //! export the type used for the permeability
    using PermeabilityType = Scalar;
    //! export the material law type used
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InjectionSpatialParams(const Problem& problem) : ParentType(problem)
    {
        layerBottom_ = 22.5;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // residual saturations
        fineMaterialParams_.setSwr(0.2);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setPe(1e4);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        return coarsePorosity_;
    }


    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    { return 790; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    { return 2700; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the solid
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    { return lambdaSolid_; }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > layerBottom_; }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    Scalar lambdaSolid_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
};

} // end namespace Dumux

#endif
