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
 * \ingroup TwoPOneCTests
 * \brief Definition of the spatial parameters for the steam injection problem
 */

#ifndef DUMUX_STEAMINJECTION_SPATIAL_PARAMS_HH
#define DUMUX_STEAMINJECTION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

//forward declaration
template<class TypeTag>
class InjectionProblemSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InjectionProblemSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InjectionProblemSpatialParams, SpatialParams, Dumux::InjectionProblemSpatialParams<TypeTag>);
}

/*!
 * \ingroup TwoPOneCTests
 * \brief Definition of the spatial parameters for various steam injection problems
 */
template<class TypeTag>
class InjectionProblemSpatialParams
: public FVSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                         typename GET_PROP_TYPE(TypeTag, Scalar),
                         InjectionProblemSpatialParams<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, InjectionProblemSpatialParams<TypeTag>>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;

public:
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = DimWorldMatrix;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InjectionProblemSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // set Van Genuchten Parameters
        materialParams_.setSwr(0.1);
        materialParams_.setSnr(0.0);
        materialParams_.setVgAlpha(0.0028);
        materialParams_.setVgn(2.0);
    }

    /*!
     * \brief Returns the hydraulic conductivity \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    DimWorldMatrix permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        DimWorldMatrix permMatrix(0.0);

        // intrinsic permeability
        permMatrix[0][0] = 1e-9;
        permMatrix[1][1] = 1e-9;

        return permMatrix; //default value
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 0.4;
    }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
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
    template<class ElementSolution>
    Scalar solidHeatCapacity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolution& elemSol) const
    { return 850.0; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 2650; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the solid
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolution& elemSol) const
    { return 2.8; }

private:
    MaterialLawParams materialParams_;
};

}

#endif
