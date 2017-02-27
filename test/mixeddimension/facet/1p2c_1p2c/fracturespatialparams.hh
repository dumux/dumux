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
#ifndef DUMUX_1P_FRACTURE_SPATIALPARAMS_HH
#define DUMUX_1P_FRACTURE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the fracture problem
 */
template<class TypeTag>
class OnePFractureSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

public:
    using PermeabilityType = Scalar;

    OnePFractureSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePermeability);
    }

    void init()
    {
        ParentType::init();

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
            return 1.0;
        else if (isBarrier(element))
            return 0.1;
        return 0.5;
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
    {
        // the barrier stores less heat (is more conductive)
        if (isBarrier(element))
            return 79;
        return 790; /*specific heat capacity of granite [J / (kg K)]*/
    }

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
    {
        // the barrier is denser
        if (isBarrier(element))
            return 3500;
        return 2700; /*density of granite [kg/m^3]*/
    }

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
    {
        // the barrier is more conductive
        if (isBarrier(element))
            return 4.5;
        return 2.8;
    }

    bool isOpenFracture(const Element& element) const
    { return isOpenFracture_[this->problem().elementMapper().index(element)]; }

    bool isBarrier(const Element& element) const
    { return isBarrier_[this->problem().elementMapper().index(element)]; }

private:

    Scalar permeability_;
    std::vector<bool> isOpenFracture_;
    std::vector<bool> isBarrier_;
};
} //end namespace

#endif
