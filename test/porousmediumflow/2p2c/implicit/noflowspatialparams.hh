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
 * \brief The spatial parameters class for the matrix problem
 */
#ifndef DUMUX_2P_NOFLOW_SPATIALPARAMS_HH
#define DUMUX_2P_NOFLOW_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag>
class TwoPNoFlowSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPSpatialParams, SpatialParams, TwoPNoFlowSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPSpatialParams, MaterialLaw)
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
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the matrix problem
 */
template<class TypeTag>
class TwoPNoFlowSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    // get the material law from the property system
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = Scalar;

    TwoPNoFlowSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        km_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixPermeability);
        kf_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePermeability);
        phim_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixPorosity);
        phif_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePorosity);
        fractureAperture_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureAperture);

        // residual saturations
        materialParamsm_.setSwr(0.0);
        materialParamsm_.setSnr(0.05);
        materialParamsf_.setSwr(0.0);
        materialParamsf_.setSnr(0.0);

        // parameters
        materialParamsm_.setPe(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixPe));
        materialParamsm_.setLambda(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixLambda));
        materialParamsf_.setPe(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePe));
        materialParamsf_.setLambda(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureLambda));

        // store info on which elements are inside the fracture
        isFractureElement_.resize(gridView.size(0), false);
        for (const auto& element : elements(gridView))
        {
            const auto p = element.geometry().center();
            const auto x = p[0];
            const auto y = p[1];
            if (y < yUpperFracture(x) && y > yLowerFracture(x))
                isFractureElement_[this->problem().elementMapper().index(element)] = true;
        }
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeability(const Element &element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        if (isFractureElement(element))
            return kf_;
        else
            return km_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosity(const Element &element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        if (isFractureElement(element))
            return phif_;
        else
            return phim_;
    }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const SubControlVolume& scv,
                                               const ElementSolutionVector& elemSol) const
    {
        if (isFractureElement(element))
            return materialParamsf_;
        else
            return materialParamsm_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    Scalar solidHeatCapacity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolutionVector& elemSol) const
    { return 790; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return 2700; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    { return 2.8; }

    bool isFractureElement(const Element& element) const
    { return isFractureElement_[this->problem().elementMapper().index(element)]; }

    Scalar yUpperFracture(Scalar x) const
    { return fractureAperture_/2; }

    Scalar yLowerFracture(Scalar x) const
    { return -fractureAperture_/2; }

private:
    Scalar km_;
    Scalar kf_;
    Scalar phim_;
    Scalar phif_;
    Scalar fractureAperture_;

    MaterialLawParams materialParamsf_;
    MaterialLawParams materialParamsm_;
    std::vector<bool> isFractureElement_;
};
} //end namespace

#endif
