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
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1PGERMANBASIN_SPATIALPARAMS_HH
#define DUMUX_1PGERMANBASIN_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePGermanBasinSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexSet = typename GridView::IndexSet;
    using ScalarVector = std::vector<Scalar>;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePGermanBasinSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView), gridView_(gridView)
    {
        int numElements = gridView_.size(0);
        paramIdx_.resize(numElements);
        porosities_.resize(numElements,0.0);
        tempInit_.resize(numElements,0.0);
        facies_.resize(numElements,0);
        modified_.resize(numElements,0);

        faciesProperties_[1] = std::tuple<Scalar,Scalar, int>(2.76, 0.0023, 1);
        faciesProperties_[2] = std::tuple<Scalar,Scalar, int>(2.5, 0.001, 2);
        faciesProperties_[10] = std::tuple<Scalar,Scalar, int>(3.57,0.003, 3);
        faciesProperties_[20] = std::tuple<Scalar,Scalar, int>(6.32, 0.004, 4);
        faciesProperties_[30] = std::tuple<Scalar,Scalar, int>(3.67, 0.00275, 5);
        faciesProperties_[40] = std::tuple<Scalar,Scalar, int>(3.57, 0.003, 6);
        faciesProperties_[41] = std::tuple<Scalar,Scalar, int>(3.67, 0.0027, 5);
        faciesProperties_[50] = std::tuple<Scalar,Scalar, int>(3.57, 0.003, 6);
        faciesProperties_[60] = std::tuple<Scalar,Scalar, int>(6.32, 0.004, 4);
        faciesProperties_[70] = std::tuple<Scalar,Scalar, int>(6.1, 0.0035, 7);
        faciesProperties_[80] = std::tuple<Scalar,Scalar, int>(3.27, 0.0012, 8);

        for (const auto& element : elements(gridView_))
        {

            //            Element eFather = element;
            //            while(eFather.level() > 0 && eFather.hasFather())
            //                eFather = eFather.father();

            int eIdx = gridView_.indexSet().index(element);
            porosities_[eIdx] = GridCreator::parameters(element)[0];
            facies_[eIdx] = GridCreator::parameters(element)[1];
            //tempInit_[eIdx] = GridCreator::parameters(element)[2];

            auto center = element.geometry().center();
            tempInit_[eIdx] =  initTemp(center);
            modified_[eIdx] = false;

            if(porosities_[eIdx] > 1.0 - 1.0e-6)
            {
                Scalar dist = 1.0e100;
                Scalar poro = 1.0;
                for(const auto& intersection : intersections(gridView_,element))
                {
                    if(intersection.neighbor())
                    {
                        const auto& elementJ = intersection.outside();
                        int eIdxJ = gridView_.indexSet().index(elementJ);
                        auto params = GridCreator::parameters(elementJ);
                        if(params[0] < 1.0 - 1.0e-6 || (porosities_[eIdxJ] < 1.0 - 1.0e-6 && porosities_[eIdxJ] > 1.0e-4))
                        {
                            auto centerJ = elementJ.geometry().center();
                            auto distVec = centerJ - center;

                            if(params[1] == facies_[eIdx])
                            {
                                auto distJ = distVec.two_norm();

                                if(distJ < dist)
                                {
                                    poro = std::max(params[0],porosities_[eIdxJ]);
                                    modified_[eIdx] = true;
                                    dist = distJ;
                                }
                            }

                        }

                    }
                }

                porosities_[eIdx] = poro;
            }

            //Scalar conductivity = GridCreator::parameters(eFather)[0];
            auto prop = faciesProperties_[facies_[eIdx]];
            //Scalar conductivity = std::pow(0.6/std::get<0>(prop),  porosities_[eIdx]) * std::get<0>(prop)/(1+std::get<1>(prop)*tempInit_[eIdx]);
            Scalar conductivity = std::pow(0.6/std::get<0>(prop),  porosities_[eIdx]) * std::get<0>(prop);
            paramIdx_[eIdx] = conductivity;

        }


        bool hasPoroOne = true;


        while(hasPoroOne)
        {
            hasPoroOne = false;
            for (const auto& element : elements(gridView_))
            {

                //            Element eFather = element;
                //            while(eFather.level() > 0 && eFather.hasFather())
                //                eFather = eFather.father();

                int eIdx = gridView_.indexSet().index(element);

                if(porosities_[eIdx] > 1.0 - 1.0e-6)
                {
                    auto center = element.geometry().center();
                    Scalar dist = 1.0e100;
                    Scalar poro = 1.0;
                    for(const auto& intersection : intersections(gridView_,element))
                    {
                        if(intersection.neighbor())
                        {
                            const auto& elementJ = intersection.outside();
                            int eIdxJ = gridView_.indexSet().index(elementJ);
                            if(porosities_[eIdxJ] < 1.0 - 1.0e-6 && porosities_[eIdxJ] > 1.0e-4)
                            {
                                auto centerJ = elementJ.geometry().center();
                                auto distVec = centerJ - center;

                                if(facies_[eIdxJ] == facies_[eIdx])
                                {
                                    auto distJ = distVec.two_norm();

                                    if(distJ < dist)
                                    {
                                        poro = porosities_[eIdxJ];
                                        modified_[eIdx] = true;
                                        dist = distJ;
                                    }
                                }

                            }

                        }
                    }

                    porosities_[eIdx] = poro;
                    if(poro > 1.0 - 1.0e-6)
                        hasPoroOne = true;
                }

                //Scalar conductivity = GridCreator::parameters(eFather)[0];
                auto prop = faciesProperties_[facies_[eIdx]];
                //Scalar conductivity = std::pow(0.6/std::get<0>(prop),  porosities_[eIdx]) * std::get<0>(prop)/(1+std::get<1>(prop)*tempInit_[eIdx]);
                Scalar conductivity = std::pow(0.6/std::get<0>(prop),  porosities_[eIdx]) * std::get<0>(prop);
                paramIdx_[eIdx] = conductivity;

            }
        }
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        int eIdx = gridView_.indexSet().index(element);

        return paramIdx_[eIdx];
    }

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    Scalar initTemp(const GlobalPosition globalPos) const
    {
        Scalar z = globalPos[dimWorld-1];
        Scalar z_max = 127.507;
        Scalar z_min = -10744.6;

        //return 281.15*(z-z_min)/(z_max-z_min) - 423.15*(z-z_max)/(z_max-z_min);
        return 281.15;
    }

public:
    std::vector<PermeabilityType> paramIdx_;
    std::vector<Scalar> porosities_;
    std::vector<Scalar> tempInit_;
    std::vector<int> facies_;
    std::vector<bool> modified_;
    //std::vector<double> paramGrid_;
    std::map<int, std::tuple<Scalar,Scalar, int>> faciesProperties_;
    const GridView gridView_;
};
} // end namespace
#endif
