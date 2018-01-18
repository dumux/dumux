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
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the column problem.
 */
#ifndef DUMUX_COLUMNXYLOL_SPATIAL_PARAMS_HH
#define DUMUX_COLUMNXYLOL_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the column problem.
 */
//forward declaration
template<class TypeTag>
class ColumnSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ColumnSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ColumnSpatialParams, SpatialParams, ColumnSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(ColumnSpatialParams, MaterialLaw)
{
 private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
 public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedParkerVanGen3P<Scalar>>;
};
} // end namespace Properties

/*!
 * \ingroup ThreePThreeCModel
 * \brief Definition of the spatial parameters for the column problem
 */
template<class TypeTag>
class ColumnSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum { dimWorld=GridView::dimensionworld };
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dimWorld>;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    ColumnSpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        // intrinsic permeabilities
        fineK_ = 1.4e-11;
        coarseK_ = 1.4e-8;

        // porosities
        finePorosity_ = 0.46;
        coarsePorosity_ = 0.46;

        // specific heat capacities
        fineHeatCap_ = 850.;
        coarseHeatCap_ = 84000.;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // residual saturations
        fineMaterialParams_.setSwr(0.12);
        fineMaterialParams_.setSnr(0.10);
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.12);
        coarseMaterialParams_.setSnr(0.10);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.5);
        fineMaterialParams_.setVgn(4.0);
        coarseMaterialParams_.setVgn(4.0);

        coarseMaterialParams_.setKrRegardsSnr(true);
        fineMaterialParams_.setKrRegardsSnr(true);

        // parameters for adsorption
        coarseMaterialParams_.setKdNAPL(0.);
        coarseMaterialParams_.setRhoBulk(1500.);
        fineMaterialParams_.setKdNAPL(0.);
        fineMaterialParams_.setRhoBulk(1500.);
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
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
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
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolutionVector& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
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
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidHeatCapacity(const Element& element,
                             const SubControlVolume& scv,
                             const ElementSolutionVector& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineHeatCap_;
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
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
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
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    {
        return lambdaSolid_;
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    {
        return (0.90 - eps_ <= globalPos[1]);
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    Scalar fineHeatCap_;
    Scalar coarseHeatCap_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    Scalar lambdaSolid_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
