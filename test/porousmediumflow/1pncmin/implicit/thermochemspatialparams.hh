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

#ifndef DUMUX_THERMOCHEM_SPATIAL_PARAMS_HH
#define DUMUX_THERMOCHEM_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/porousmediumflow/1pncmin/implicit/indices.hh>
#include <dumux/material/fluidmatrixinteractions/porosityreactivebed.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ThermoChemSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ThermoChemSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ThermoChemSpatialParams, SpatialParams, Dumux::ThermoChemSpatialParams<TypeTag>);

} // end namespace Properties

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxTestProblems
 * \brief Definition of the spatial parameters for the FuelCell
 *        problem which uses the isothermal 2p2c box model
 */
template<class TypeTag>
class ThermoChemSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ThisType = ThermoChemSpatialParams<TypeTag>;
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using Element = typename GridView::template Codim<0>::Entity;

    using PorosityLaw = PorosityReactiveBed<TypeTag>;
    using PermeabilityLaw = PermeabilityKozenyCarman<TypeTag>;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    ThermoChemSpatialParams(const Problem& problem, const GridView &gridView)
    : ParentType(problem, gridView)
    {
        //thermal conductivity of CaO
        lambdaSolid_ = 0.4; //[W/(m*K)] Nagel et al [2013b]

        eps_ = 1e-6;
    }

    ~ThermoChemSpatialParams()
    {}

    /*!
     * \brief Called by the Problem to initialize the spatial params.
     */
    void init()
    {
        //! Intitialize the parameter laws
        poroLaw_.init(*this);
        permLaw_.init(*this);
    }

    /*!
     *  \brief Define the initial permeability \f$[m^2]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar initialPermeability(const Element& element, const SubControlVolume &scv) const
    { return 5e-12; }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
//     Scalar porosity(const Element &element,
//                     const FVElementGeometry &fvGeometry,
//                     const int scvIdx) const
//     {
//         const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
//
//         if (globalPos[1]<eps_)
//             return porosity_;
//         else
//             return 0.2;
//
//     }


    /*!
     *  \brief Define the initial porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar initialPorosity(const Element& element, const SubControlVolume &scv) const
    {
         Scalar phi;

         if(isCharge_==true)
         phi = 0.8;  //direct charging acc. to Nagel et al 2014
//              phi = 0.887; //indirect charging acc. to Schmitt 2016
         else
          phi = 0.604;  //direct charging acc. to Nagel et al 2014
//             phi = 0.772;  //indirect charging acc. to Schmitt 2016

         return phi;
    }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *
     *  Solution dependent permeability function
     */
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return permLaw_.evaluatePermeability(element, scv, elemSol); }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minPorosity(const Element& element, const SubControlVolume &scv) const
    { return 0.604; //intrinsic porosity of CaO2H2 (see Nagel et al. 2014)
        //        return 0.772;  //indirect charging acc. to Schmitt 2016
    }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ after clogging caused by mineralization
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    { return poroLaw_.evaluatePorosity(element, scv, elemSol); }



    Scalar solidity(const SubControlVolume &scv) const
    { return 1.0 - porosityAtPos(scv.center()); }

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
        return 790;
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
//      return 3370; // density of CaO [kg/m^3]
        return 2600;
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
    { return lambdaSolid_; }

private:

   Scalar eps_;
   Scalar lambdaSolid_;
   bool isCharge_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, IsCharge);
   PorosityLaw poroLaw_;
   PermeabilityLaw permLaw_;
};

}//end namespace

#endif
