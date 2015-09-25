// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_LENSSPATIALPARAMETERS_2P_HH
#define DUMUX_LENSSPATIALPARAMETERS_2P_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/boxmodels/2p/2pmodel.hh>

/**
 * @file
 * @brief Class for defining spatial parameters
 * @author Bernd Flemisch, Klaus Mosthaf, Markus Wolff
 */

namespace Dumux
{

/** \todo Please doc me! */

template<class TypeTag>
class LensSpatialParameters2p : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<CoordScalar,dimWorld,dimWorld> FieldMatrix;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    // define the material law which is parameterized by effective
    // saturations
//    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    LensSpatialParameters2p(const GridView& gridView)
        : ParentType(gridView),
          lensK_(0),
          outerK_(0)
    {
        Dumux::InterfaceSoilProperties interfaceSoilProps("interface2p.xml");

        lensPorosity_ = interfaceSoilProps.ISP_FinePorosity;
        outerPorosity_ = interfaceSoilProps.ISP_CoarsePorosity;

        for(int i = 0; i < dim; i++){
            lensK_[i][i] = interfaceSoilProps.ISP_FinePermeability;
            outerK_[i][i] = interfaceSoilProps.ISP_CoarsePermeability;
        }

        // residual saturations
        lensMaterialParams_.setSwr(interfaceSoilProps.ISP_FineResidualSaturationWetting);
        lensMaterialParams_.setSnr(interfaceSoilProps.ISP_FineResidualSaturationNonWetting);
        outerMaterialParams_.setSwr(interfaceSoilProps.ISP_CoarseResidualSaturationWetting);
        outerMaterialParams_.setSnr(interfaceSoilProps.ISP_CoarseResidualSaturationNonWetting);

        // parameters for the Van Genuchten law
        // alpha and n
//        lensMaterialParams_.setVgAlpha(0.00045);
//        lensMaterialParams_.setVgN(7.3);
//        outerMaterialParams_.setVgAlpha(0.0037);
//        outerMaterialParams_.setVgN(4.7);

        // parameters for the Brooks-Corey law
        lensMaterialParams_.setPe(interfaceSoilProps.ISP_FineBrooksCoreyEntryPressure);
        lensMaterialParams_.setAlpha(interfaceSoilProps.ISP_FineBrooksCoreyLambda);
        outerMaterialParams_.setPe(interfaceSoilProps.ISP_CoarseBrooksCoreyEntryPressure);
        outerMaterialParams_.setAlpha(interfaceSoilProps.ISP_CoarseBrooksCoreyLambda);

        // parameters for the linear law
        // minimum and maximum pressures
 //        lensMaterialParams_.setEntryPC(0);
//        outerMaterialParams_.setEntryPC(0);
//        lensMaterialParams_.setMaxPC(0);
//        outerMaterialParams_.setMaxPC(0);
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    const FieldMatrix& intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isInLens_(globalPos))
            return lensPorosity_;
        return outerPorosity_;
    }

    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;

        if (isInLens_(globalPos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }


    //! Set the bounding box of the fine-sand lens
    void setLensCoords(const GlobalPosition& lensLowerLeft,
                       const GlobalPosition& lensUpperRight)
    {
        lensLowerLeft_ = lensLowerLeft;
        lensUpperRight_ = lensUpperRight;
    }

private:
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    FieldMatrix lensK_;
    FieldMatrix outerK_;
    Scalar lensPorosity_;
    Scalar outerPorosity_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;
};

} // end namespace
#endif

