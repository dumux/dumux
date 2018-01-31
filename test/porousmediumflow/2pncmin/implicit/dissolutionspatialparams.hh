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
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
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

/*!
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
template<class TypeTag>
class DissolutionSpatialparams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLawParams = typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using CoordScalar = typename GridView::ctype;
    enum { dimWorld=GridView::dimensionworld };

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using PorosityLaw = PorosityPrecipitation<TypeTag>;
    using PermeabilityLaw = PermeabilityKozenyCarman<TypeTag>;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Tensor;

    DissolutionSpatialparams(const Problem& problem)
    : ParentType(problem)
    {
        solubilityLimit_     = getParam<Scalar>("SpatialParams.SolubilityLimit", 0.26);
        referencePorosity_     = getParam<Scalar>("SpatialParams.referencePorosity", 0.11);
        refPermeability_ = getParam<Scalar>("SpatialParams.referencePermeability", 2.23e-14);
        irreducibleLiqSat_   = getParam<Scalar>("SpatialParams.IrreducibleLiqSat", 0.2);
        irreducibleGasSat_   = getParam<Scalar>("SpatialParams.IrreducibleGasSat", 1e-3);
        pEntry1_             = getParam<Scalar>("SpatialParams.Pentry1", 500);
        bcLambda1_           = getParam<Scalar>("SpatialParams.BCLambda1", 2);

        // residual saturations
        materialParams_.setSwr(irreducibleLiqSat_);
        materialParams_.setSnr(irreducibleGasSat_);

        // parameters of Brooks & Corey Law
        materialParams_.setPe(pEntry1_);
        materialParams_.setLambda(bcLambda1_);

        //! Intitialize the parameter laws
        poroLaw_.init(*this);
        permLaw_.init(*this);

        for (int i = 0; i < dimWorld; i++)
            referencePermeability_[i][i] = refPermeability_;
    }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minPorosity(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }

    /*!
     *  \brief Define the reference porosity \f$[-]\f$ distribution.
     *  This is the porosity of the porous medium without any of the
     *  considered solid phases.
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return referencePorosity_; }

    /*!
     *  \brief Return the actual recent porosity \f$[-]\f$ accounting for
     *  clogging caused by mineralization
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    { return poroLaw_.evaluatePorosity(element, scv, elemSol); }

    /*!
     *  \brief Define the reference permeability \f$[m^2]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    PermeabilityType referencePermeability(const Element& element, const SubControlVolume &scv) const
    { return referencePermeability_; }


    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *
     *  Solution dependent permeability function
     */
    PermeabilityType permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return permLaw_.evaluatePermeability(element, scv, elemSol); }

    Scalar solidity(const SubControlVolume &scv) const
    { return 1.0 - porosityAtPos(scv.center()); }

    Scalar solubilityLimit() const
    { return solubilityLimit_; }

    Scalar theta(const SubControlVolume &scv) const
    { return 10.0; }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

private:

    MaterialLawParams materialParams_;

    PorosityLaw poroLaw_;
    PermeabilityLaw permLaw_;
    Scalar solubilityLimit_;
    Scalar referencePorosity_;
    Scalar refPermeability_;
    PermeabilityType referencePermeability_ = 0.0;
    Scalar irreducibleLiqSat_;
    Scalar irreducibleGasSat_;
    Scalar pEntry1_;
    Scalar bcLambda1_;
};

} // end namespace Dumux

#endif
