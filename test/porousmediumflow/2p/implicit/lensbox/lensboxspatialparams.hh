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
 * \brief The spatial parameters for the SatConditionProblem which uses the
 *        two-phase fully implicit model
 */
#ifndef DUMUX_LENS_SPATIAL_PARAMS_HH
#define DUMUX_LENS_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/2p/implicit/vertextominpcelemmapper.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class SatConditionSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(SatConditionSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(SatConditionSpatialParams, SpatialParams, Dumux::SatConditionSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(SatConditionSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    //vangenuchten
    //typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
    // brookscorey
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;//RawMaterialLaw;
public:
    // define the material law parameterized by absolute saturations
    //vangenuchten
    //typedef EffToAbsLaw<EffectiveLaw> type;
    //brookscorey
    typedef EffToAbsLaw<EffectiveLaw> type;//EffToAbsLaw<RawMaterialLaw> type;
};
}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the SatConditionProblem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class SatConditionSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef typename Dumux::VertIdxToMinPcMapper<TypeTag> VertIdxToMinPcMapper;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    SatConditionSpatialParams(const GridView& gridView)
    : ParentType(gridView)
    {
            lensLowerLeft_[0]   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensLowerLeftX);
            lensLowerLeft_[1]   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensLowerLeftY);
            lensUpperRight_[0]  = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensUpperRightX);
            lensUpperRight_[1]  = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensUpperRightY);

        // residual saturations
        fineMaterialParams_.setSwr(0.098);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.078);
        coarseMaterialParams_.setSnr(0.0);

        //// parameters of Brooks & Corey Law
        fineMaterialParams_.setPe(1324);
        fineMaterialParams_.setLambda(2.49);
        coarseMaterialParams_.setPe(370);
        coarseMaterialParams_.setLambda(3.86);

        // parameters for the Van Genuchten law
        // alpha and n
        //fineMaterialParams_. setVgAlpha(0.000581);
        //fineMaterialParams_.setVgn(5.34);
        //coarseMaterialParams_.setVgAlpha(0.00225);
        //coarseMaterialParams_.setVgn(8.06);

        fineK_ = 5.26e-11;
        coarseK_ = 5.04e-10;
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {
        const GlobalPosition& globalPos = element.geometry().center();

        if (isFine_(globalPos))
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
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    { const GlobalPosition& globalPos = element.geometry().center();
        if (isFine_(globalPos))
            return 0.35;
        return 0.4;
        }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
    {
        const GlobalPosition& globalPos = element.geometry().center();
        if (isFine_(globalPos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

private:
    bool isFine_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] || globalPos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    Scalar fineK_;
    Scalar coarseK_;
    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
};

} // end namespace
#endif
