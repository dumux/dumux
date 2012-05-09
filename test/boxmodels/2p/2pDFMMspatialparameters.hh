// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
/*!
 * \file
 *
 * \brief The spatial parameters for the LensProblem which uses the
 *        twophase box model
 */
#ifndef DUMUX_2PDFMM_SPATIAL_PARAMETERS_HH
#define DUMUX_2PDFMM_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/boxmodels/2pDFMM/2pdfmmodel.hh>

#include "dumux/io/artreader_new.hh"

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPDFMMSpatialParameters;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPDFMMSpatialParameters);

// Set the spatial parameters
SET_TYPE_PROP(TwoPDFMMSpatialParameters, SpatialParams, Dumux::TwoPDFMMSpatialParameters<TypeTag>);

// Set the material Law
SET_PROP(TwoPDFMMSpatialParameters, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}
/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxTestProblems
 * \brief The spatial parameters for the LensProblem which uses the
 *        twophase box model
 */
template<class TypeTag>
class TwoPDFMMSpatialParameters : public BoxSpatialParameters<TypeTag>
{

	template<int dim>
	struct FaceLayout
	{
		bool contains (Dune::GeometryType gt)
		{
			return gt.dim() == dim - 1;
		}
	};
	template<int dim>
	struct VertexLayout
	{
		bool contains (Dune::GeometryType gt)
		{
			return gt.dim() == 0;
		}
	};



    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,VertexLayout> VM;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,FaceLayout> FM;



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

    TwoPDFMMSpatialParameters(const GridView& gridView)
        : ParentType(gridView)
    {
        try
        {
//            lensLowerLeft_[0]   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensLowerLeftX);
//            lensLowerLeft_[1]   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensLowerLeftY);
//            lensUpperRight_[0]  = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensUpperRightX);
//            lensUpperRight_[1]  = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensUpperRightY);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

        int nEdges = fractureEdgesIdx_.size();
        vFractureWidth.resize(nEdges);

//        // residual saturations
//        lensMaterialParams_.setSwr(0.18);
//        lensMaterialParams_.setSnr(0.0);
//        outerMaterialParams_.setSwr(0.05);
//        outerMaterialParams_.setSnr(0.0);
//
//        // parameters for the Van Genuchten law
//        // alpha and n
//        lensMaterialParams_.setVgAlpha(0.00045);
//        lensMaterialParams_.setVgN(7.3);
//        outerMaterialParams_.setVgAlpha(0.0037);
//        outerMaterialParams_.setVgN(4.7);




//         parameters for the linear law
//         minimum and maximum pressures
//        lensMaterialParams_.setEntryPC(0);
//        outerMaterialParams_.setEntryPC(0);
//        lensMaterialParams_.setMaxPC(0);
//        outerMaterialParams_.setMaxPC(0);

        lensK_  		= 9.05e-12;
        outerK_ 		= 4.6e-10;

    	gridView_ 		= 0;
    	facemapper_ 	= 0;
    	vertexmapper_ 	= 0;
    	setupFractureMatrixSoilParameters();

    }



    void setupFractureMatrixSoilParameters()
    {
    	Scalar mD = 1e-12 * 1e-3; //miliDarcy

    	SwrF_ 	 = 0.00;
    	SwrM_ 	 = 0.00;
    	SnrF_ 	 = 0.00;
    	SnrM_ 	 = 0.00;
    	pdF_  	 = 2.5*1e4;
    	pdM_   	 = 2.5*1e4;
    	lambdaF_ = 2.0;
    	lambdaM_ = 2.0;

    	rockMatrixMaterialParams_.setSwr(SwrM_);
        rockMatrixMaterialParams_.setSnr(SnrM_);
        fractureMaterialParams_.setSwr(SwrF_);
        fractureMaterialParams_.setSnr(SnrF_);

        rockMatrixMaterialParams_.setPe(pdM_);
        rockMatrixMaterialParams_.setLambda(lambdaM_);
        fractureMaterialParams_.setPe(pdF_);
        fractureMaterialParams_.setLambda(lambdaF_);

        KMatrix_   = 1 * mD; //m^2
        KFracture_ = 1e4 * mD; //8.3 * 1e7 * mD; //m^2

        porosityMatrix_   =	0.25;
        porosityFracture_ = 0.10;
        fractureWidth_    = 1e-1;

        for (int i=0; i<vFractureWidth.size();i++ )
        {
        	vFractureWidth[i]  = 1e-1;
        }
    }



    /*!
     * \brief Intrinsic permability
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Intrinsic permeability
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    Scalar intrinsicPermeabilityFracture(const Element &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int scvIdx) const
    {
    	int globalIdx = vertexmapper().map(element, scvIdx, dim);
        if (isDuneFractureVertex_[globalIdx])
            return KFracture_;
        else
        return KFracture_;
//        DUNE_THROW(Dune::InvalidStateException,
//                "Called soil::KFracture() for non-fracture vertex");
    }
    /*!
     * \brief Porosity
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Porosity
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    { return porosityMatrix_; }

    // TODO correlate with function fractureWidth
    Scalar porosityFracture(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
//    	return 0.5;//homogeneous
    	return porosityFracture_;
    }
    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;

        if (isInLens_(globalPos))
            return rockMatrixMaterialParams_;
        return rockMatrixMaterialParams_;
    }

    const MaterialLawParams& materialLawParamsFracture(const Element &element,
                                                    const FVElementGeometry &fvElemGeom,
                                                    int scvIdx) const
	{
//        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
		int globalIdx = vertexmapper().map(element, scvIdx, dim);

		// be picky if called for non-fracture vertices
		assert(isVertexFracture(globalIdx));

		return fractureMaterialParams_;
	}


    bool isVertexFracture(const Element &element, int localVertexIdx) const
    {
        if (inactivate_fractures == true)
        	return false;
        else {
    	int globalIdx = vertexmapper().map(element, localVertexIdx, dim);
        return isDuneFractureVertex_[globalIdx];
        }
    }

    bool isVertexFracture(int globalIdx) const
    {
        if (inactivate_fractures == true)
        	return false;
        else {
        	return isDuneFractureVertex_[globalIdx];
        }
    }

    bool isEdgeFracture(const Element &element, int localFaceIdx) const
    {
        int globalIdx = facemapper().map(element, localFaceIdx, 1);
        return isDuneFractureEdge_[globalIdx];
    }

    const VM &vertexmapper() const
    { return *vertexmapper_;};

    const FM &facemapper() const
    { return *facemapper_;};

    Scalar fractureWidth(int globalFaceIdx) const
    {
    	return fractureWidth_;
//    	return vFractureWidth(globalFaceIdx); //TODO to be implemented
    }

    Scalar fractureWidth(const Element &element, int localFaceIdx) const
    {
    	int globalFaceIdx = facemapper().map(element, localFaceIdx, 1);
//    	return vFractureWidth[globalFaceIdx]; //TODO to be implemented
    	return fractureWidth_;
    }


    void setGridView(const GridView &gv)
     {
       delete gridView_;
       delete facemapper_;
       delete vertexmapper_;

       gridView_ = new GridView(gv);
       facemapper_ = new FM(gv);
       vertexmapper_ = new VM(gv);
     }


    const void setFractureBoolVectors(std::vector<bool>& isDuneFractureVertex,
                                std::vector<bool>& isDuneFractureEdge,
                                std::vector<int>& fractureEdgesIdx,
                                bool useFractures)
    {
        isDuneFractureVertex_ 	= isDuneFractureVertex;
        isDuneFractureEdge_ 	= isDuneFractureEdge;
        inactivate_fractures 	= useFractures;
        fractureEdgesIdx_    	= fractureEdgesIdx;
/*
        std::cout << "isDuneFractureVertex_ =";
        for (int i = 0; i < isDuneFractureVertex_.size(); ++i) {
            std::cout << isDuneFractureVertex_[i] << "  ";
        };
        std::cout << "\n";

        std::cout << "isDuneFractureEdge_ =";
        for (int i = 0; i < isDuneFractureEdge_.size(); ++i) {
            std::cout << isDuneFractureEdge_[i] << "  ";
        };
        std::cout << "\n";
*/
    }

    Scalar SwrF_, SwrM_, SnrF_, SnrM_, lambdaF_, lambdaM_, pdF_, pdM_;

private:
    bool isInLens_(const GlobalPosition &pos) const
    {
//        for (int i = 0; i < dim; ++i) {
//            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
//                return false;
//        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;

    Scalar KMatrix_;
    Scalar KFracture_;
    Scalar porosityMatrix_;
    Scalar porosityFracture_;

    std::vector<int> vFractureWidth;
    Scalar fractureWidth_;

    MaterialLawParams fractureMaterialParams_;
    MaterialLawParams rockMatrixMaterialParams_;
    bool inactivate_fractures;

    std::vector<bool> isFracture_;
    std::vector<bool> isDuneFractureVertex_;
    std::vector<bool> isDuneFractureEdge_;
    std::vector<int>  fractureEdgesIdx_;

    VM *vertexmapper_;
    FM *facemapper_;
    GridView *gridView_;


};

} // end namespace
#endif

