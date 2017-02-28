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
 * \brief The spatial parameters class for the fracture problem
 */
#ifndef DUMUX_2P2C_FRACTURE_SPATIALPARAMS_HH
#define DUMUX_2P2C_FRACTURE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPTwoCFractureSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPTwoCFractureSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPTwoCFractureSpatialParams, SpatialParams, TwoPTwoCFractureSpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(TwoPTwoCFractureSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}

template<class TypeTag>
class TwoPTwoCFractureSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

public:
    using PermeabilityType = Scalar;

    TwoPTwoCFractureSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePermeability);

    }

    void init()
    {
        ParentType::init();

        // residual saturations
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.0);
        materialLawParamsBarrier_.setSwr(0.2);
        materialLawParamsBarrier_.setSnr(0.0);
        materialLawParamsConduit_.setSwr(0.2);
        materialLawParamsConduit_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        materialLawParams_.setPe(1e4);
        materialLawParams_.setLambda(2.0);
        materialLawParamsBarrier_.setPe(1e5);
        materialLawParamsBarrier_.setLambda(2.0);
        materialLawParamsConduit_.setPe(1e2);
        materialLawParamsConduit_.setLambda(2.0);

        // initialize which elements are open fractures/barriers
        const auto numElements = this->problem().gridView().size(0);
        isOpenFracture_.resize(numElements, false);
        isBarrier_.resize(numElements, false);
        for (const auto& element : elements(this->problem().gridView()))
        {
            const auto globalPos = element.geometry().center();
            if (this->problem().couplingManager().bulkProblem().isOpenFracture(globalPos))
                isOpenFracture_[this->problem().elementMapper().index(element)] = true;
            else if (this->problem().couplingManager().bulkProblem().isBarrier(globalPos))
                isBarrier_[this->problem().elementMapper().index(element)] = true;
        }
    }

    /*!
     * \brief Return the intrinsic permeability for a given position in [m^2].
     */
    Scalar permeability(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        using std::pow;
        static const Scalar openK = pow(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureAperture), 2)/12;
        static const Scalar barrierK = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixPermeability)/50;

        if (isOpenFracture(element))
            return openK;
        else if (isBarrier(element))
            return barrierK;
        return permeability_;
    }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The solution for all dofs of the element
     */
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return 0; }

    /*!
     * \brief Define the porosity in [-].
     */
    Scalar porosity(const Element &element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        if (isOpenFracture(element))
            return 0.8;
        else if (isBarrier(element))
            return 0.1;
        return 0.5;
    }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const SubControlVolume& scv,
                                               const ElementSolutionVector& elemSol) const
    {
        if (isOpenFracture(element))
            return materialLawParamsConduit_;
        else if (isBarrier(element))
            return materialLawParamsBarrier_;
        return materialLawParams_;
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

    bool isOpenFracture(const Element& element) const
    { return isOpenFracture_[this->problem().elementMapper().index(element)]; }

    bool isBarrier(const Element& element) const
    { return isBarrier_[this->problem().elementMapper().index(element)]; }

private:

    Scalar permeability_;
    MaterialLawParams materialLawParams_;
    MaterialLawParams materialLawParamsBarrier_;
    MaterialLawParams materialLawParamsConduit_;

    std::vector<bool> isOpenFracture_;
    std::vector<bool> isBarrier_;
};
} //end namespace

#endif
