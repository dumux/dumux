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
 * \ingroup OnePTests
 * \brief Definition of the spatial parameters for the 1pni problems.
 */
#ifndef DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH
#define DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTests
 * \brief Definition of the spatial parameters for the 1pni problems.
 */
template<class TypeTag>
class OnePNISpatialParams
: public FVSpatialParamsOneP<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                             typename GET_PROP_TYPE(TypeTag, Scalar),
                             OnePNISpatialParams<TypeTag>>
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, OnePNISpatialParams<TypeTag>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePNISpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-10; }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    { return 790; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    { return 2700; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    { return 2.8; }
};

} // end namespace Dumux

#endif
