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
 * \brief Definition of the spatial parameters for the kuevette problem.
 */
#ifndef DUMUX_KUEVETTE3P3CNI_SPATIAL_PARAMS_HH
#define DUMUX_KUEVETTE3P3CNI_SPATIAL_PARAMS_HH

#include <dune/common/float_cmp.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the kuevette problem.
 */
//forward declaration
template<class TypeTag>
class KuevetteSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(KuevetteSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(KuevetteSpatialParams, SpatialParams, KuevetteSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(KuevetteSpatialParams, MaterialLaw)
{
 private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
 public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedParkerVanGen3P<Scalar>>;
};
} // end namespace Dumux

/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the kuevette problem
 */
template<class TypeTag>
class KuevetteSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum { dimWorld=GridView::dimensionworld };

    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
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
    KuevetteSpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        // intrinsic permeabilities
        fineK_ = 6.28e-12;
        coarseK_ = 9.14e-10;

        // porosities
        finePorosity_ = 0.42;
        coarsePorosity_ = 0.42;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // residual saturations
        fineMaterialParams_.setSwr(0.12);
        fineMaterialParams_.setSnr(0.07);
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.12);
        coarseMaterialParams_.setSnr(0.07);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.005);
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
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    {
        return 850; // specific heat capacity of sand [J / (kg K)]
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
        return ((Dune::FloatCmp::ge<Scalar>(globalPos[0], 0.13)
                    && Dune::FloatCmp::le<Scalar>(globalPos[0], 1.24)
                    && Dune::FloatCmp::ge<Scalar>(globalPos[1], 0.32)
                    && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.60))
                || (Dune::FloatCmp::ge<Scalar>(globalPos[0], 1.20)
                    && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.15)));
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    Scalar lambdaSolid_;
};

} // end namespace Dumux

#endif
