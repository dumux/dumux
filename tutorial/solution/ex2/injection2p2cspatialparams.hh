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
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
// TODO: dumux-course-task
// Inlcude your own material law
#include "mymateriallaw.hh"

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotmateriallaw.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag>
class InjectionSpatialParams;

// setup property TypeTag
namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InjectionSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InjectionSpatialParams, SpatialParams, InjectionSpatialParams<TypeTag>);

// TODO: dumux-course-task
// Use your own material law instead
// Set the material law parameterized by absolute saturations
SET_PROP(InjectionSpatialParams, MaterialLaw)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    // using type = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
    using type = EffToAbsLaw<MyMaterialLaw<Scalar>>;
};

} // end namespace Properties

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the injection problem
 *        which uses the isothermal two-phase two-component
 *        fully implicit model.
 */
template<class TypeTag>
class InjectionSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param problem The problem
     */
    InjectionSpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        aquiferHeightFromBottom_ = 30.0;

        // intrinsic permeabilities
        aquitardK_ = getParam<Scalar>("SpatialParams.PermeabilityAquitard");
        aquiferK_ = getParam<Scalar>("SpatialParams.PermeabilityAquifer");

        // porosities
        aquitardPorosity_ = 0.2;
        aquiferPorosity_ = 0.4;

        // residual saturations
        aquitardMaterialParams_.setSwr(0.2);
        aquitardMaterialParams_.setSnr(0.0);
        aquiferMaterialParams_.setSwr(0.2);
        aquiferMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        aquitardMaterialParams_.setPe(getParam<Scalar>("SpatialParams.EntryPressureAquitard"));
        aquiferMaterialParams_.setPe(getParam<Scalar>("SpatialParams.EntryPressureAquifer"));
        aquitardMaterialParams_.setLambda(2.0);
        aquiferMaterialParams_.setLambda(2.0);

        // plot the material laws using gnuplot and exit
        if (getParam<bool>("Problem.OnlyPlotMaterialLaws"))
        {
            plotMaterialLaws();
            exit(0);
        }
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const

    {
        if (isInAquitard_(globalPos))
            return aquitardK_;
        return aquiferK_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInAquitard_(globalPos))
            return aquitardPorosity_;
        return aquiferPorosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param globalPos The global position
     *
     * \return the material parameters object
     */
     const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isInAquitard_(globalPos))
            return aquitardMaterialParams_;
        return aquiferMaterialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

    /*!
     * \brief Creates a gnuplot output of the pc-Sw curve
     */
    void plotMaterialLaws()
    {
        PlotMaterialLaw<Scalar, MaterialLaw> plotMaterialLaw;
        GnuplotInterface<Scalar> gnuplot;
        plotMaterialLaw.addpcswcurve(gnuplot, aquitardMaterialParams_, 0.2, 1.0, "upper layer (fine, aquitard)", "w lp");
        plotMaterialLaw.addpcswcurve(gnuplot, aquiferMaterialParams_, 0.2, 1.0, "lower layer (coarse, aquifer)", "w l");
        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("pc-Sw");
    }

private:

    static constexpr Scalar eps_ = 1e-6;

    bool isInAquitard_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > aquiferHeightFromBottom_ + eps_; }

    Scalar aquitardK_;
    Scalar aquiferK_;
    Scalar aquiferHeightFromBottom_;

    Scalar aquitardPorosity_;
    Scalar aquiferPorosity_;

    MaterialLawParams aquitardMaterialParams_;
    MaterialLawParams aquiferMaterialParams_;
};

} // end namespace Dumux

#endif
