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
#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/2pncmin/implicit/indices.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class DissolutionSpatialparams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(DissolutionSpatialparams);

// Set the spatial parameters
SET_TYPE_PROP(DissolutionSpatialparams, SpatialParams, DissolutionSpatialparams<TypeTag>);

// Set the material Law
SET_PROP(DissolutionSpatialparams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
};

} // end namespace Properties

/**
 * \brief Definition of the spatial parameters for the brine-co2 problem
 *
 */
template<class TypeTag>
class DissolutionSpatialparams : public ImplicitSpatialParams<TypeTag>
{
    using ThisType = DissolutionSpatialparams<TypeTag>;
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLawParams = typename GET_PROP_TYPE(TypeTag, MaterialLawParams);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using CoordScalar = typename GridView::ctype;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
    };

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using Element = typename GridView::template Codim<0>::Entity;

    using PorosityLaw = PorosityPrecipitation<TypeTag>;
    using PermeabilityLaw = PermeabilityKozenyCarman<TypeTag, Scalar>;

public:
    DissolutionSpatialparams(const Problem& problem, const GridView &gridView)
    : ParentType(problem, gridView),
      permLaw_(problem)
    {
        // residual saturations
        materialParams_.setSwr(0.2);
        materialParams_.setSnr(1e-3);

        // parameters of Brooks & Corey Law
        materialParams_.setPe(500);
        materialParams_.setLambda(2);
    }

    /*!
     * \brief Called by the Problem to initialize the spatial params.
     */
    void init()
    {
        using namespace std::placeholders;
        using InitParamFunction = std::function<Scalar(const Element&, const SubControlVolume&)>;

        InitParamFunction minPoro = std::bind(&ThisType::porosityMin, this, _1, _2);
        InitParamFunction initPoro = [](const Element&, const SubControlVolume&) { return 0.11; };
        InitParamFunction initPerm = [](const Element&, const SubControlVolume&) { return 2.23e-14; };

        poroLaw_.init(std::move(initPoro), std::move(minPoro));
        initPoro = [](const Element&, const SubControlVolume&) { return 0.11; };
        permLaw_.init(std::move(initPoro), std::move(initPerm));
    }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Solution dependent permeability function
     */
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return permLaw_.evaluatePermeability(element, scv, elemSol); }

    /*!
     * \brief Define the minimum porosity \f$[-]\f$ after salt precipitation
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosityMin(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }

    /*!
     * \brief Define the minimum porosity \f$[-]\f$ after clogging caused by mineralization
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    { return poroLaw_.evaluatePorosity(element, scv, elemSol); }


    Scalar solidity(const SubControlVolume &scv) const
    { return 1.0 - porosityAtPos(scv.center()); }

    Scalar solubilityLimit() const
    { return 0.26; }

    Scalar theta(const SubControlVolume &scv) const
    { return 10.0; }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

private:
    MaterialLawParams materialParams_;

    PorosityLaw poroLaw_;
    PermeabilityLaw permLaw_;
};

} // end namespace Dumux

#endif
