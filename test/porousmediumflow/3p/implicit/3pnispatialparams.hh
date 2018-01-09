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
 * \ingroup ThreePTests
 * \brief Definition of the spatial parameters for the 3pni problems.
 */
#ifndef DUMUX_THREEPNI_SPATIAL_PARAMS_HH
#define DUMUX_THREEPNI_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePTests
 * \brief Definition of the spatial parameters for the 3pni problems.
 */

//forward declaration
template<class TypeTag>
class ThreePNISpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ThreePNISpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ThreePNISpatialParams, SpatialParams, ThreePNISpatialParams<TypeTag>);

// Set the material Law
SET_PROP(ThreePNISpatialParams, MaterialLaw)
{
 private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = RegularizedParkerVanGen3P<Scalar>;
 public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};
}


template<class TypeTag>
class ThreePNISpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum { dimWorld=GridView::dimensionworld };
    using CoordScalar = typename Grid::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    ThreePNISpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSnr(0.10);
        materialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(0.5);
        materialParams_.setVgn(4.0);
        materialParams_.setKrRegardsSnr(true);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    {
        return 790; // specific heat capacity of granite [J / (kg K)]
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    {
        return 2700; // density of granite [kg/m^3]
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    {
        return lambdaSolid_;
    }


private:

    MaterialLawParams materialParams_;
    Scalar permeability_;
    Scalar porosity_;
    Scalar lambdaSolid_;
};

} // end namespace Dumux

#endif
