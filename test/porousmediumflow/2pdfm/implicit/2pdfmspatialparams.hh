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
 * \brief The spatial parameters for the 2pDFM Problem which uses the
 *        twophase discrete fracture model.
 */
#ifndef DUMUX_TEST_2PDFM_SPATIAL_PARAMETERS_HH
#define DUMUX_TEST_2PDFM_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/2pdfm/implicit/model.hh>
#include <dumux/io/artgridcreator.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <fstream>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPDFMSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPDFMSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPDFMSpatialParams, SpatialParams, Dumux::TwoPDFMSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPDFMSpatialParams, MaterialLaw)
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
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the 2PDFMProblem which uses the
 *        twophase model
 */
template<class TypeTag>
class TwoPDFMSpatialParams : public ImplicitSpatialParams<TypeTag>
{

    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType geomType)
        {
            return geomType.dim() == dim - 1;
        }
    };

    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, FaceLayout> FaceMapper;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    TwoPDFMSpatialParams(const GridView& gridView)
        : ParentType(gridView), gridView_(gridView),
        faceMapper_(gridView), vertexMapper_(gridView),
        fractureMapper_(gridView)
    {
        useSimpleTest_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, UseSimpleTest);
        inactivateFractures_ = false;

        Scalar mD = 1e-12 * 1e-3; //miliDarcy

        rockMatrixMaterialParams_.setSwr(0.1);
        rockMatrixMaterialParams_.setSnr(0.0);

        fineFractureMaterialParams_.setSwr(0.09);
        fineFractureMaterialParams_.setSnr(0.0);
        coarseFractureMaterialParams_.setSwr(0.08);
        coarseFractureMaterialParams_.setSnr(0.0);

        rockMatrixMaterialParams_.setPe(2000);
        rockMatrixMaterialParams_.setLambda(2.49);

        fineFractureMaterialParams_.setPe(1300);
        fineFractureMaterialParams_.setLambda(2.0);
        coarseFractureMaterialParams_.setPe(300);
        coarseFractureMaterialParams_.setLambda(3.2);

        KMatrix_   = 1 * mD; //m^2
        KFracture_ = 1e5 * mD; //m^2

        porosityMatrix_   = 0.25;
        porosityFracture_ = 0.10;
        fractureWidth_    = 1e-2;

        // comment this out if you want to use the simple test case
        fractureMapper_.map();

        isVertexFracture_.resize(gridView.size(dim));
        for (const auto& element : elements(gridView))
        {
            for (int scvIdx = 0; scvIdx < element.subEntities(dim); scvIdx++)
            {
                int vIdxGlobal = vertexMapper_.subIndex(element, scvIdx, dim);
                isVertexFracture_[vIdxGlobal] = isVertexFracture(element, scvIdx);
            }
        }
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
                                 int scvIdx) const
    {
        return KMatrix_;
    }

    /*!
     * \brief Intrinsic permeability of fractures.
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     */
    Scalar intrinsicPermeabilityFracture(const Element &element,
                                         const FVElementGeometry &fvGeometry,
                                         int scvIdx) const
    {
        return KFracture_;
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
                    int scvIdx) const
    { return porosityMatrix_; }

    /*!
     * \brief Porosity Fracture
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Porosity Fracture
     */
    Scalar porosityFracture(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        return porosityFracture_;
    }
    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
    {
        return rockMatrixMaterialParams_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object of the Fracture
     */
    const MaterialLawParams& materialLawParamsFracture(const Element &element,
                                                    const FVElementGeometry &fvGeometry,
                                                    int scvIdx) const
    {
        DUNE_UNUSED int vIdxGlobal = vertexMapper_.subIndex(element, scvIdx, dim);

        // be picky if called for non-fracture vertices
        assert(isVertexFracture(vIdxGlobal));

        const GlobalPosition& globalPos = element.geometry().center();
        if (isFine_(globalPos))
            return fineFractureMaterialParams_;
        return coarseFractureMaterialParams_;
    }

    /*!
     * \brief Checks whether vertex is a fracture.
     *
     * \param element The current element
     * \param localVertexIdx Vertex index to be checked
     */
    bool isVertexFracture(const Element &element, int localVertexIdx) const
    {
        if (inactivateFractures_)
        {
            return false;
        }

        if (!useSimpleTest_)
        {
            int vIdxGlobal = vertexMapper_.subIndex(element, localVertexIdx, dim);
            return fractureMapper_.isDuneFractureVertex(vIdxGlobal);
        }

        const auto vertexPosition = element.template subEntity<dim>(localVertexIdx).geometry().center();
        return (vertexPosition[1] > 0.5 - 1e-6 && vertexPosition[1] < 0.5 + 1e-6) ||
               (vertexPosition[0] > 0.5 - 1e-6 && vertexPosition[0] < 0.5 + 1e-6 && vertexPosition[1] > 0.25 && vertexPosition[1] < 0.75);
    }

    /*!
     * \brief Checks whether vertex is a fracture.
     *
     * \param vIdxGlobal Vertex index to be checked
     */
    bool isVertexFracture(int vIdxGlobal) const
    {
        if (inactivateFractures_)
        {
            return false;
        }

        if (!useSimpleTest_)
            return fractureMapper_.isDuneFractureVertex(vIdxGlobal);
        return isVertexFracture_[vIdxGlobal];
    }

    /*!
     * \brief Checks whether element edge is a fracture.
     *
     * \param element The current element
     * \param localFaceIdx Face index to be checked
     */
    bool isEdgeFracture(const Element &element, int localFaceIdx) const
    {
        if (!useSimpleTest_)
        {
            int fIdxGlobal = faceMapper_.subIndex(element, localFaceIdx, 1);
            return fractureMapper_.isDuneFractureEdge(fIdxGlobal);
        }
        const auto edgePosition = element.template subEntity<1>(localFaceIdx).geometry().center();
        return (edgePosition[1] > 0.5 - 1e-6 && edgePosition[1] < 0.5 + 1e-6) ||
               (edgePosition[0] > 0.5 - 1e-6 && edgePosition[0] < 0.5 + 1e-6 && edgePosition[1] > 0.25 && edgePosition[1] < 0.75);
    }

    bool hasFractureFaces(const Element& element, const FVElementGeometry& fvGeometry, int localVertexIdx) const
    {
        // determine whether or not the scv has fracture faces connected to it
        for (int fIdx = 0; fIdx < element.subEntities(1); ++fIdx)
        {
          auto& scvFace = fvGeometry.subContVolFace[fIdx];
          if (isEdgeFracture(element, fIdx) && (scvFace.i == localVertexIdx || scvFace.j == localVertexIdx))
              return true;
        }
        return false;
    }

    /*!
     * \brief Returns the width of the fracture.
     *
     * \param globalFaceIdx Global face index of which the width is returned
     */
    Scalar fractureWidth(int globalFaceIdx) const
    {
        return fractureWidth_;
    }

    /*!
     * \brief Returns the width of the fracture.
     *
     * \param element The current element
     * \param localFaceIdx Local face index of which the width is returned
     */
    Scalar fractureWidth(const Element &element, int localFaceIdx) const
    {
        return fractureWidth_;
    }

private:
    bool isFine_(const GlobalPosition &globalPos) const
    {
           if (globalPos[0] > 0.5)
               return true;
           else
               return false;
    }

    bool useSimpleTest_;

    Scalar KMatrix_;
    Scalar KFracture_;
    Scalar porosityMatrix_;
    Scalar porosityFracture_;

    Scalar fractureWidth_;

    MaterialLawParams coarseFractureMaterialParams_;
    MaterialLawParams fineFractureMaterialParams_;
    MaterialLawParams rockMatrixMaterialParams_;
    bool inactivateFractures_;

    const GridView gridView_;
    const FaceMapper faceMapper_;
    const VertexMapper vertexMapper_;

    std::vector<bool> isVertexFracture_;
    Dumux::FractureMapper<TypeTag> fractureMapper_;
};

} // end namespace
#endif // DUMUX_TEST_2PDFM_SPATIAL_PARAMETERS_HH
