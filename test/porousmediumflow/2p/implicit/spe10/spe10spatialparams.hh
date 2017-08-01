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
#ifndef DUMUX_SPE10_SPATIAL_PARAMS_HH
#define DUMUX_SPE10_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include "quadraticlaw.hh"

#include "spe10permeability.hh"
#include "spe10porosity.hh"

namespace Dumux
{

//forward declaration
template<class TypeTag>
class Spe10SpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Spe10SpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Spe10SpatialParams, SpatialParams, Dumux::Spe10SpatialParams<TypeTag>);

// Set the material Law
SET_PROP(Spe10SpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef QuadraticLaw<Scalar,RegularizedBrooksCorey<Scalar>> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};

}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the SPE10 test problem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class Spe10SpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dimWorld,dimWorld> DimWorldMatrix;

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
    Spe10SpatialParams(const GridView& gridView)
    : ParentType(gridView)
    {
        if (dimWorld == 2)
        {
                layerIdx_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, LayerIdx);
                omittedCoordinate_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, OmittedCoordinate);
        }

        auto numElements = gridView.size(0);

        porosity_.resize(numElements);
        permeability_.resize(numElements);
        materialLawParams_.resize(numElements);

        const unsigned num3DElements = 1122000;

        Scalar oneMDarcyInM2 = 9.869233e-16;

        DimWorldMatrix meanPermeability(0);
        Scalar meanPorosity = 0.0;
        Scalar totalVolume = 0.0;

        for (const auto& element : elements(gridView))
        {
            unsigned int idx = getIndex_(element.geometry().center());
            unsigned int idxGrid = GridCreator::grid().leafGridView().indexSet().index(element);

            permeability_[idxGrid] = 0;
            for (int i = 0; i < dimWorld; i++)
                permeability_[idxGrid][i][i] = oneMDarcyInM2*SPE10Permeability<Scalar>::data[i*num3DElements + idx];

            porosity_[idxGrid] = std::max(SPE10Porosity<Scalar>::data[idx],1.0e-12);

            Scalar volume = element.geometry().volume();
            totalVolume += volume;

            for (int i = 0; i < dim; i++)
            {
                meanPermeability[i][i] += volume / permeability_[idxGrid][i][i];
            }
            meanPorosity += volume * porosity_[idxGrid];
        }


        for (int i = 0; i < dim; i++)
        {
            meanPermeability[i][i] /= totalVolume;
            meanPermeability[i][i] = 1.0 / meanPermeability[i][i];
        }

        meanPorosity /= totalVolume;

        Scalar referencePd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, ReferenceEntryPressure);

        for (int eIdx = 0; eIdx < numElements; eIdx++)
        {
            Scalar elementEntryPressure = referencePd * sqrt((porosity_[eIdx]*meanPermeability[1][1]) / (meanPorosity * permeability_[eIdx][1][1]));
            if (porosity_[eIdx] < 1e-8)
                elementEntryPressure = referencePd * sqrt((meanPermeability[1][1]) / (meanPorosity * permeability_[eIdx][1][1]));
            materialLawParams_[eIdx].setPe(elementEntryPressure);
            materialLawParams_[eIdx].setLambda(2.0);
            materialLawParams_[eIdx].setSwr(0.0);
            materialLawParams_[eIdx].setSnr(0.0);
        }
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const DimWorldMatrix intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return permeability_[eIdx];
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
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return porosity_[eIdx];
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
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return materialLawParams_[eIdx];
    }

private:
    unsigned int getIndex_(const GlobalPosition& globalPos) const
    {
        int i,j,k;
        Scalar hX = 6.096;
        Scalar hY = 3.048;
        Scalar hZ = 0.6096;

        if (dim == 2)
        {
            switch (omittedCoordinate_)
            {
              case 0:
                i = layerIdx_;
                j = globalPos[0]/hY;
                k = globalPos[1]/hZ;
                break;
              case 1:
                i = globalPos[0]/hX;
                j = layerIdx_;
                k = globalPos[1]/hZ;
                break;
              default:
                i = globalPos[0]/hX;
                j = globalPos[1]/hY;
                k = layerIdx_;
                break;
            }
        }
        else
        {
            i = globalPos[0]/hX;
            j = globalPos[1]/hY;
            k = globalPos[2]/hZ;
        }

        int nX = 60;
        int nY = 220;
        return k*(nX*nY) + j*nX + i;
    }

    std::vector<DimWorldMatrix> permeability_;
    std::vector<Scalar> porosity_;
    std::vector<MaterialLawParams> materialLawParams_;
    int layerIdx_;
    int omittedCoordinate_;
};

} // end namespace
#endif
