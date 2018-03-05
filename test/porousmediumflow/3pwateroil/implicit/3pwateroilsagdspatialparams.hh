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
 * \ingroup ThreePWaterOilTests
 * \brief Definition of the spatial parameters for the SAGD problem.
 */

#ifndef DUMUX_SAGD_SPATIAL_PARAMS_HH
#define DUMUX_SAGD_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>

#include <dumux/porousmediumflow/3pwateroil/indices.hh>

#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilTests
 * \brief Definition of the spatial parameters for the SAGD problem.
 */
//forward declaration
template<class TypeTag>
class SagdSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(SagdSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(SagdSpatialParams, SpatialParams, Dumux::SagdSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(SagdSpatialParams, MaterialLaw)
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

/*!
 * \ingroup ThreePThreeCNIModel
 *
 * \brief Definition of the spatial parameters for the SAGD problem
 */
template<class TypeTag>
class SagdSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    enum { dimWorld=GridView::dimensionworld };

    using CoordScalar = typename Grid::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    SagdSpatialParams(const Problem& problem)
    : ParentType(problem), eps_(1e-6)
    {
        layerBottom_ = 35.0;

        // intrinsic permeabilities
        fineK_ = 1e-16;
        coarseK_ = 4e-14;

        // porosities
        finePorosity_ = 0.10;
        coarsePorosity_ = 0.1;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // specific heat capacities
        fineHeatCap_ = 850.;
        coarseHeatCap_ = 850;

        // residual saturations
        fineMaterialParams_.setSwr(0.1);
        fineMaterialParams_.setSwrx(0.12);  //Total liquid Residual Saturation
        fineMaterialParams_.setSnr(0.09);   //Residual of NAPL if there is no water
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.1);
        coarseMaterialParams_.setSwrx(0.12);
        coarseMaterialParams_.setSnr(0.09);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgn(4.0);
        coarseMaterialParams_.setVgn(4.0);
        fineMaterialParams_.setVgAlpha(1.);
        coarseMaterialParams_.setVgAlpha(1.);

        coarseMaterialParams_.setKrRegardsSnr(false);
        fineMaterialParams_.setKrRegardsSnr(false);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {  return permeabilityAtPos(scv.dofPosition());}

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
        else
            return coarsePorosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the material parameters object
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        return materialLawParamsAtPos(scv.dofPosition());
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
        else
            return coarseMaterialParams_;
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
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineHeatCap_ ;
        else
            return coarseHeatCap_;
    }

     /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    {
        return 2650; // density of sand [kg/m^3]
    }

     /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    {
        return lambdaSolid_;
    }


private:
    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        return pos[dimWorld-1] > layerBottom_ - eps_;
    };

    Scalar layerBottom_;
    Scalar lambdaSolid_;

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    Scalar fineHeatCap_;
    Scalar coarseHeatCap_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    Scalar eps_;
};

}

#endif
