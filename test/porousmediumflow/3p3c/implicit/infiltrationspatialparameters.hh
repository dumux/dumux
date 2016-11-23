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
 * \brief Definition of the spatial parameters for the kuevette problem
 *        which uses the isothermal two-phase two component
 *        fully implicit model.
 */
#ifndef DUMUX_INFILTRATION_SPATIAL_PARAMETERS_HH
#define DUMUX_INFILTRATION_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/3p3c/implicit/indices.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class InfiltrationSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InfiltrationSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InfiltrationSpatialParams, SpatialParams, InfiltrationSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(InfiltrationSpatialParams, MaterialLaw)
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
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class TypeTag>
class InfiltrationSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim = GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> DimVector;


    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
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
    InfiltrationSpatialParams(const Problem& problem, const GridView &gridView)
        : ParentType(problem, gridView)
    {
        // intrinsic permeabilities
        fineK_ = 1.e-11;
        coarseK_ = 1.e-11;

        // porosities
        porosity_ = 0.40;

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(0.0005);
        materialParams_.setVgn(4.);
        materialParams_.setKrRegardsSnr(false);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar intrinsicPermeability(const SubControlVolume &scv,
                                 const VolumeVariables& volVars) const
    {
        if (isFineMaterial_(scv.dofPosition()))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }


    /*!
     * \brief return the parameter object for the material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return
            70.0 <= globalPos[0] && globalPos[0] <= 85.0 &&
            7.0 <= globalPos[1] && globalPos[1] <= 7.50;
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;
};

}

#endif
