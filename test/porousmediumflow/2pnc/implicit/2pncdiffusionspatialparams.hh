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
 * \ingroup TwoPNCTest
 * \brief Definition of the spatial parameters for the fuel cell
 *        problem which uses the isothermal/non-insothermal 2pnc box model
 */

#ifndef DUMUX_TWOPNC_DIFFUSION_SPATIAL_PARAMS_HH
#define DUMUX_TWOPNC_DIFFUSION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPNCDiffusionSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPNCDiffusionSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPNCDiffusionSpatialParams, SpatialParams, TwoPNCDiffusionSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPNCDiffusionSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};

} // end namespace Properties

/*!
 * \ingroup TwoPNCTest
 * \brief Definition of the spatial parameters for the TwoPNCDiffusion
 *        problem which uses the isothermal 2p2c box model
 */
template<class TypeTag>
class TwoPNCDiffusionSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using CoordScalar = typename GridView::ctype;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = DimWorldMatrix;

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    TwoPNCDiffusionSpatialParams(const Problem& problem)
    : ParentType(problem), K_(0)
    {
        // intrinsic permeabilities
        K_[0][0] = 5e-11;
        K_[1][1] = 5e-11;

        // porosities
        porosity_ = 0.2;

        // residual saturations
        materialParams_.setSwr(0.02); //here water, see philtophoblaw
        materialParams_.setSnr(0.0);

        //parameters for the vanGenuchten law
        materialParams_.setPe(1e2); // alpha = 1/pcb
        materialParams_.setLambda(1.3);
    }

    /*!
     * \brief Returns the hydraulic conductivity \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    DimWorldMatrix permeabilityAtPos(const GlobalPosition& globalPos) const
    { return K_; }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:
    DimWorldMatrix K_;
    Scalar porosity_;
    static constexpr Scalar eps_ = 1e-6;
    MaterialLawParams materialParams_;
};

}//end namespace

#endif
