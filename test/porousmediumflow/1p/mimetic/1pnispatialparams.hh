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
 * \brief Definition of the spatial parameters for the 1pni problems.
 */
#ifndef DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH
#define DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

#include "../../2p/mimetic/spe10/spe10permeability.hh"
#include "../../2p/mimetic/spe10/spe10porosity.hh"


namespace Dumux
{

/*!
 * \ingroup OnePNIModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the 1pni problems.
 */
template<class TypeTag>
class OnePNISpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
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

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    OnePNISpatialParams(const Problem& problem, const GridView &gridView)
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
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    { return 790; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    { return 2700; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    { return 2.8; }

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
    int layerIdx_;
    int omittedCoordinate_;

};

} // end namespace Dumux

#endif
