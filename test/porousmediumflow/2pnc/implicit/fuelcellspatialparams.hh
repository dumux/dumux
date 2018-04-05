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
 * \ingroup TwoPNCTests
 * \brief Definition of the spatial parameters for the fuel cell
 *        problem which uses the isothermal/non-insothermal 2pnc box model
 */

#ifndef DUMUX_FUELCELL_SPATIAL_PARAMS_HH
#define DUMUX_FUELCELL_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/philtophoblaw.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCTests
 * \brief Definition of the spatial parameters for the fuel cell
 *        problem which uses the isothermal/non-insothermal 2pnc box model
 */
//forward declaration
template<class TypeTag>
class FuelCellSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(FuelCellSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(FuelCellSpatialParams, SpatialParams, FuelCellSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(FuelCellSpatialParams, MaterialLaw)
{
 private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffMaterialLaw = RegularizedVanGenuchten<Scalar>;

 public:
    // define the material law parameterized by absolute saturations
    using type = PhilToPhobLaw<EffMaterialLaw>;
};

} // end namespace Properties

/*!
 * \ingroup TwoPNCMinModel
 * \ingroup BoxTestProblems
 * \brief Definition of the spatial parameters for the FuelCell
 *        problem which uses the isothermal 2p2c box model
 */
template<class TypeTag>
class FuelCellSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using CoordScalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = DimWorldMatrix;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    FuelCellSpatialParams(const Problem& problem)
    : ParentType(problem), K_(0)
    {
        // intrinsic permeabilities
        K_[0][0] = 5e-11;
        K_[1][1] = 5e-11;

        // porosities
        porosity_ = 0.2;

        // residual saturations
        materialParams_.setSwr(0.12); //here water, see philtophoblaw
        materialParams_.setSnr(0.0);

        //parameters for the vanGenuchten law
        materialParams_.setVgAlpha(6.66e-5); // alpha = 1/pcb
        materialParams_.setVgn(3.652);
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
    {
        if (globalPos[1] < eps_)
            return porosity_;
        else
            return 0.2;
    }

    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

private:
    DimWorldMatrix K_;
    Scalar porosity_;
    static constexpr Scalar eps_ = 1e-6;
    MaterialLawParams materialParams_;
};

}//end namespace

#endif
