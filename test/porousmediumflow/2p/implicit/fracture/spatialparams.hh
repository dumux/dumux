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
 * \ingroup TwoPTests
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
#ifndef DUMUX_TWOP_FRACTURE_TEST_SPATIALPARAMS_HH
#define DUMUX_TWOP_FRACTURE_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class FractureSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(FractureSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(FractureSpatialParams, SpatialParams, Dumux::FractureSpatialParams<TypeTag>);

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
 * \ingroup TwoPTests
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class FractureSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    enum { dimWorld=GridView::dimensionworld };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! export permeability type
    using PermeabilityType = Scalar;
    //! export the type used for the material law
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    FractureSpatialParams(const Problem& problem): ParentType(problem)
    {
        // residual saturations
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(0.0037);
        materialParams_.setVgn(4.7);
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-10; }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The position at which we evaluate
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    MaterialLawParams materialParams_;
};

} // end namespace Dumux

#endif
