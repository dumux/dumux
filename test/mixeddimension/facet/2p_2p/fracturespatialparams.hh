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
#ifndef DUMUX_2P_FRACTURE_SPATIALPARAMS_HH
#define DUMUX_2P_FRACTURE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag>
class TwoPFractureSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(FractureSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(FractureSpatialParams, SpatialParams, TwoPFractureSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(FractureSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the fracture problem
 */
template<class TypeTag>
class TwoPFractureSpatialParams : public ImplicitSpatialParams<TypeTag>
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

public:
    using PermeabilityType = Scalar;

    TwoPFractureSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePermeability);

        // residual saturations
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(0.0037);
        materialParams_.setVgn(4.7);
    }

    /*!
     * \brief Return the intrinsic permeability for a given position in [m^2].
     */
    Scalar permeability(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return permeability_; }

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
    { return 0.8; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

private:

    Scalar permeability_;
    MaterialLawParams materialParams_;
};
} //end namespace

#endif
