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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief The spatial parameters for the LensProblem which uses the
 *        twophase box model
 */
#ifndef DUMUX_TWOPMINC_SPATIAL_PARAMS_HH
#define DUMUX_TWOPMINC_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>

#include <dumux/implicit/2pminc/2pmincmodel.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPMincSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPMincSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPMincSpatialParams, SpatialParams, Dumux::TwoPMincSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPMincSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief The spatial parameters for the 2pMinc Test Problem which uses the
 *        twophase box model
 */
template<class TypeTag>
class TwoPMincSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum {
        numEq       = GET_PROP_VALUE(TypeTag,NumEq),
        numContinua = GET_PROP_VALUE(TypeTag, NumContinua),

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

    TwoPMincSpatialParams(const GridView& gridView)
        : ParentType(gridView)
    {
        fractureMaterialParams_.setSwr(0.00);
        fractureMaterialParams_.setSnr(0.00);
        matrixMaterialParams_.setSwr(0.00);
        matrixMaterialParams_.setSnr(0.00);

        //parameters for the Brooks-Corey law
        fractureMaterialParams_.setPe(1000.0);
        fractureMaterialParams_.setLambda(2.0);
        matrixMaterialParams_.setPe(2000.0);
        matrixMaterialParams_.setLambda(2.0);

        KFracture_ = 1.0e-8;

        for (int nC = 0; nC<numContinua; nC++)
        {
            if (nC == 0)
                KMatrix_[nC] = KFracture_;
            else
                KMatrix_[nC] = 1.0E-12; //homogeneous
        }

        porosityFracture_   = 0.9;
        porosityMatrix_     = 0.3;
    }

    /*!
     * \brief Intrinsic permability
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Intrinsic permeability
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx, int nC) const
    {
        if (nC == 0)
            return KFracture_;
        
        return KMatrix_[nC];
    }

    /*!
     * \brief Porosity
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Porosity
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx, int nC) const
    {
        if (nC == 0)
            return porosityFracture_;
        
        return porosityMatrix_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParamsFracture(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
    {
        return fractureMaterialParams_;
    }
    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsMatrix(const Element &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int scvIdx) const
    {
        return matrixMaterialParams_;
    }

private:
    Scalar KMatrix_[numContinua];
    Scalar KFracture_;

    MaterialLawParams fractureMaterialParams_;
    MaterialLawParams matrixMaterialParams_;
    Scalar porosityFracture_;
    Scalar porosityMatrix_;
};

} // end namespace
#endif
