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
#ifndef DUMUX_MIMETIC_SPE10_SPATIAL_PARAMS_HH
#define DUMUX_MIMETIC_SPE10_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
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
SET_TYPE_PROP(Spe10SpatialParams, SpatialParams, Spe10SpatialParams<TypeTag>);

// Set the material Law
SET_PROP(Spe10SpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = QuadraticLaw<Scalar,RegularizedBrooksCorey<Scalar>>;
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};
}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class Spe10SpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    //get the material law from the property system
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    Spe10SpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
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
                permeability_[idxGrid][i][i] = oneMDarcyInM2*(SPE10Permeability<Scalar>::data[i*num3DElements + idx]);

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
     * \param globalPos The global position
     */
    PermeabilityType  permeability(const Element& element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return permeability_[eIdx];
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return porosity_[eIdx];
    }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    {
        return materialLawParams_[0];
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

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
