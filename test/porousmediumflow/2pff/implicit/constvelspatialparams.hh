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
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
#ifndef DUMUX_FRACFLOW_CONSTVEL_SPATIAL_PARAMS_HH
#define DUMUX_FRACFLOW_CONSTVEL_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedtestlaw.hh>
#include <dumux/io/plotmateriallaw.hh>
#include<dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ConstVelSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ConstVelSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ConstVelSpatialParams, SpatialParams, ConstVelSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(ConstVelSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = Regularizedtestlaw<Scalar>;
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};
}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the ConstVelProblem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class ConstVelSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    //get the material law from the property system
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    ConstVelSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        // residual saturations
       materialParams_.setSwr(0.0);
       materialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        // alpha and n
        materialParams_.setPe(5171.068);
        materialParams_.setLambda(4.0);

        K_ = 100 * 9.86923* 1.0e-16;

PlotMaterialLaw<TypeTag> plotMaterialLaw;
        GnuplotInterface<Scalar> gnuplot(true);

      plotMaterialLaw.addDswcurve( gnuplot, materialParams_, 0.0, 1.0);
       gnuplot.setOption("set xrange [0:1]");
       gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
       gnuplot.plot("D-Sw");



   /*    plotMaterialLaw.addlamswcurve( gnuplot, materialParams_, 0.0, 1.0);

        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("LambdarTerm-Sw");




plotMaterialLaw.addkrcurves( gnuplot, materialParams_, 0.0, 1.0);

        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("Kr-Sw");689


plotMaterialLaw.addpcswcurve( gnuplot, materialParams_, 0.0, 1.0);

        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("pc-Sw");


       plotMaterialLaw.adddpcdswcurve( gnuplot, materialParams_, 0.0, 1.0);
       gnuplot.setOption("set xrange [0.1:1]");
       gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
       gnuplot.plot("dpcdsw-Sw");
*/

    }

 /*!

     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return K_;
    }

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
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
    }

private:
    Scalar K_;
    MaterialLawParams materialParams_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
