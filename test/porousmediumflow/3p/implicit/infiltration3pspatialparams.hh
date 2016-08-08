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
 * \brief Definition of the spatial parameters for the kuevette problem, which
 *        uses the three-phase fully implicit model.
 */
#ifndef DUMUX_INFILTRATION_THREEP_SPATIAL_PARAMS_HH
#define DUMUX_INFILTRATION_THREEP_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/3p3c/implicit/indices.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>
#include <dumux/io/plotmateriallaw3p.hh>
namespace Dumux
{

//forward declaration
template<class TypeTag>
class InfiltrationThreePSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InfiltrationThreePSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InfiltrationThreePSpatialParams, SpatialParams, InfiltrationThreePSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(InfiltrationThreePSpatialParams, MaterialLaw)
{
 private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedParkerVanGen3P<Scalar> EffectiveLaw;
 public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup ThreePModel
 *
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class TypeTag>
class InfiltrationThreePSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;


public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InfiltrationThreePSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        // intrinsic permeabilities
        fineK_ = GET_RUNTIME_PARAM(TypeTag, Scalar, permeability);
        coarseK_ = GET_RUNTIME_PARAM(TypeTag, Scalar, permeability);

        // porosities
        porosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, porosity);

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(GET_RUNTIME_PARAM(TypeTag, Scalar, vanGenuchtenAlpha));
        materialParams_.setVgn(GET_RUNTIME_PARAM(TypeTag, Scalar, vanGenuchtenN));
        materialParams_.setKrRegardsSnr(false);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);

        plotFluidMatrixInteractions_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output,
                                                                    PlotFluidMatrixInteractions);
    }

    ~InfiltrationThreePSpatialParams()
    {}

     /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        PlotMaterialLaw<TypeTag> plotMaterialLaw(plotFluidMatrixInteractions_);

        plotMaterialLaw.plotpc(materialParams_);
        plotMaterialLaw.plotkr(materialParams_);
    }

    /*!
     * \brief Intrinsic permability
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     * \return Intrinsic permeability
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Porosity
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     * \return Porosity
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return porosity_;
    }


    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvElemGeom,
                                               int scvIdx) const
    {
        return materialParams_;
    }

    const MaterialLawParams& materialLawParams() const
    {
        return materialParams_;
    }
private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return
            70. <= globalPos[0] && globalPos[0] <= 85. &&
            7.0 <= globalPos[1] && globalPos[1] <= 7.50;
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;

    bool plotFluidMatrixInteractions_;
};

}

#endif
