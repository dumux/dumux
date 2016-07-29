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
 * \brief Definition of the spatial parameters for the fuel cell
 *        problem which uses the isothermal/non-insothermal 2pnc box model
 */

#ifndef DUMUX_FUELCELL_SPATIAL_PARAMS_HH
#define DUMUX_FUELCELL_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/philtophoblaw.hh>
#include <dumux/porousmediumflow/2pnc/implicit/model.hh>

namespace Dumux
{

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
class FuelCellSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename Grid::ctype;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,

        wPhaseIdx = FluidSystem::wPhaseIdx
    };

    using GlobalPosition = Dune::FieldVector<CoordScalar,dimWorld>;
    using DimVector = Dune::FieldVector<CoordScalar,dim>;
    using DimMatrix = Dune::FieldMatrix<CoordScalar,dim,dim>;
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    FuelCellSpatialParams(const Problem& problem, const GridView &gridView)
    : ParentType(problem, gridView), K_(0)
    {
        // intrinsic permeabilities
        K_[0][0] = 5e-11;
        K_[1][1] = 5e-11;

        // porosities
        porosity_ = 0.2;

        //thermalconductivity
        lambdaSolid_ = 14.7; //[W/(m*K)] Acosta et al. [2006]

        // residual saturations
        materialParams_.setSwr(0.12); //here water, see philtophoblaw
        materialParams_.setSnr(0.0);

        //parameters for the vanGenuchten law
        materialParams_.setVgAlpha(6.66e-5); // alpha = 1/pcb
        materialParams_.setVgn(3.652);

        eps_ = 1e-6;
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const DimMatrix intrinsicPermeability(const SubControlVolume& scv) const
    { return K_; }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosity(const SubControlVolume& scv) const
    {
        const auto& globalPos = scv.dofPosition();

        if (globalPos[1]<eps_)
            return porosity_;
        else
            return 0.2;

    }

    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const SubControlVolume& scv) const
    {
        return materialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar heatCapacity(const Element &element,
                        const SubControlVolume& scv) const
    {
        return
            790 // specific heat capacity of granite [J / (kg K)]
            * 2700 // density of granite [kg/m^3]
            * (1 - porosity(scv));
    }

    // /*!
    //  * \brief Calculate the heat flux \f$[W/m^2]\f$ through the
    //  *        rock matrix based on the temperature gradient \f$[K / m]\f$
    //  *
    //  * This is only required for non-isothermal models.
    //  *
    //  * \param heatFlux The resulting heat flux vector
    //  * \param fluxVars The flux variables
    //  * \param elemVolVars The volume variables
    //  * \param tempGrad The temperature gradient
    //  * \param element The current finite element
    //  * \param fvGeometry The finite volume geometry of the current element
    //  * \param faceIdx The local index of the sub-control volume face where
    //  *                    the matrix heat flux should be calculated
    //  */
    // void matrixHeatFlux(DimVector &heatFlux,
    //                     const FluxVariables &fluxVars,
    //                     const ElementVolumeVariables &elemVolVars,
    //                     const DimVector &tempGrad,
    //                     const Element &element,
    //                     const FVElementGeometry &fvGeometry,
    //                     const int faceIdx) const
    // {

    //     static const Scalar lWater = 0.6;
    //     static const Scalar lGranite = 2.8;

    //     // arithmetic mean of the liquid saturation and the porosity
    //     const int i = fvGeometry.subContVolFace[faceIdx].i;
    //     const int j = fvGeometry.subContVolFace[faceIdx].j;
    //     Scalar sW = std::max<Scalar>(0.0, (elemVolVars[i].saturation(wPhaseIdx) +
    //                                        elemVolVars[j].saturation(wPhaseIdx)) / 2);
    //     Scalar poro = (porosity(element, fvGeometry, i) +
    //                    porosity(element, fvGeometry, j)) / 2;

    //     Scalar lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
    //     Scalar ldry = pow(lGranite, (1-poro));

    //     // the heat conductivity of the matrix. in general this is a
    //     // tensorial value, but we assume isotropic heat conductivity.
    //     Scalar heatCond = ldry + sqrt(sW) * (ldry - lsat);

    //     // the matrix heat flux is the negative temperature gradient
    //     // times the heat conductivity.
    //     heatFlux = tempGrad;
    //     heatFlux *= -heatCond;
    // }

    Scalar thermalConductivitySolid(const Element &element,
                                    const SubControlVolume& scv) const
    {
        return lambdaSolid_;
    }

private:
    DimMatrix K_;
    Scalar porosity_;
    Scalar eps_;
    MaterialLawParams materialParams_;
    Scalar lambdaSolid_;
};

}//end namespace

#endif
