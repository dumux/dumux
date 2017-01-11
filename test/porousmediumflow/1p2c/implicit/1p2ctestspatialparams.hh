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
 * \brief Definition of the spatial parameters for the 1p2c
 *        outlfow problem.
 */
#ifndef DUMUX_1P2C_TEST_SPATIAL_PARAMS_HH
#define DUMUX_1P2C_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the 1p2c
 *        outflow problem.
 */
template<class TypeTag>
class OnePTwoCTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePTwoCTestSpatialParams(const Problem& problem, const GridView &gridView)
        : ParentType(problem, gridView)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

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
    { return lambdaSolid_; }



private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar lambdaSolid_;
};

} // end namespace Dumux

#endif
