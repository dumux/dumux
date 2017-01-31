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
 * \brief Definition of the spatial parameters for the effective diffusivity tests.
 */
#ifndef DUMUX_FLUIDMATRIXINTERACTION_TEST_SPATIAL_PARAMS_HH
#define DUMUX_FLUIDMATRIXINTERACTION_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class FluidMatrixInteractionTestSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(FluidMatrixInteractionTestSpatialParams);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(FluidMatrixInteractionTestSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Define whether to open a gnuplot window
NEW_PROP_TAG(OutputOpenPlotWindow);
SET_BOOL_PROP(FluidMatrixInteractionTestSpatialParams, OutputOpenPlotWindow, false);
}

/*!
 * \ingroup MaterialTestProblems
 * \brief Definition of the spatial parameters for the effective diffusivity tests.
 */
template<class TypeTag>
class FluidMatrixInteractionTestSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:

    //! exort permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    FluidMatrixInteractionTestSpatialParams(const Problem& problem, const GridView &gridView)
        : ParentType(problem, gridView)
    {
        porosity_ = 0.3;
        rhoSolid_ = 2700.0;
        lambdaSolid_ = 2.8;

        // residual saturations
        materialParams_.setSwr(0.2);
        materialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        materialParams_.setPe(1e4);
        materialParams_.setLambda(2.0);
    }

    //! \copydoc ImplicitProblem::permeabilityAtPos()
    PermeabilityType permeabilityAtPos(const GlobalPosition &globalPos) const
    { return 1e-10; }

    //! \copydoc ImplicitProblem::porosityAtPos()
    Scalar porosityAtPos(const GlobalPosition &globalPos) const
    { return porosity_; }

    //! \copydoc ImplicitProblem::materialLawParamsAtPos()
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    { return materialParams_; }

    //! \copydoc ImplicitProblem::solidHeatCapacityAtPos()
    Scalar solidHeatCapacityAtPos(const GlobalPosition &globalPos) const
    { return 790; /* specific heat capacity of granite [J / (kg K)] */ }

    //! \copydoc ImplicitProblem::solidDensityAtPos()
    Scalar solidDensityAtPos(const GlobalPosition &globalPos) const
    { return rhoSolid_; /* density of granite [kg/m^3] */ }

    //! \copydoc ImplicitProblem::solidThermalConductivityAtPos()
    Scalar solidThermalConductivityAtPos(const GlobalPosition &globalPos) const
    { return lambdaSolid_; /* [W/(m K) */ }

protected:
    Scalar porosity_;
    Scalar lambdaSolid_;
    Scalar rhoSolid_;

    MaterialLawParams materialParams_;
};

} // end namespace Dumux

#endif
