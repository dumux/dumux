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
 * \brief The spatial parameters class for the root system test problem
 */
#ifndef DUMUX_ROOTSYSTEM_TEST_SPATIALPARAMS_HH
#define DUMUX_ROOTSYSTEM_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the root system test problem
 */
template<class TypeTag>
class RootsystemTestSpatialParams: public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld
    };
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    struct RootParams
    {
        Scalar radius;
        Scalar surface;
        Scalar axialPerm;
        Scalar radialPerm;
        int order;
        int branchId;
        Scalar mass;
    };

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootsystemTestSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView), gridView_(gridView)
    {
        Kx_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Kx);
        Kr_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Kr);

        rootParams_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.gridView().indexSet().index(element);
            auto level0element = element;
            for(auto levelIdx = element.level(); levelIdx != 0; levelIdx--)
                level0element = level0element.father();
            Scalar rootLength = element.geometry().volume();
            Scalar rootSurface = GridCreator::parameters(level0element)[2]/(1 << element.level());

            rootParams_[eIdx].radius = rootSurface / rootLength / 2.0 / M_PI;
            rootParams_[eIdx].order = GridCreator::parameters(level0element)[0];
            // root branch  -> count from 0!!
            rootParams_[eIdx].branchId = GridCreator::parameters(level0element)[1] -1;
            rootParams_[eIdx].surface = rootSurface;
            rootParams_[eIdx].mass = GridCreator::parameters(level0element)[3];

            if ((int)rootParams_[eIdx].order == 1)
            {
                rootParams_[eIdx].axialPerm = Kx_; //Kx
                rootParams_[eIdx].radialPerm = Kr_; //Kr
            }
            else if  ((int)rootParams_[eIdx].order == 2)
            {
                rootParams_[eIdx].axialPerm = Kx_;  //Kx
                rootParams_[eIdx].radialPerm = Kr_;  //Kr
            }
            else //order >= 3
                rootParams_[eIdx].axialPerm = Kx_; //Kx
                rootParams_[eIdx].radialPerm = Kr_; //Kr
        }
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param ipGlobal The integration point
     * \note Kx has units [m^4] so we have to divide by the cross-section area
     */
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        auto eIdx = gridView_.indexSet().index(element);
        const auto r = radius(eIdx);
        return Kx_ / (M_PI*r*r);
    }

    /*!
     * \brief Return the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param the index of the element
     */
    Scalar radius(unsigned int eIdxGlobal) const
    {
        return rootParams_[eIdxGlobal].radius;
    }

    Scalar rootSurface(unsigned int eIdxGlobal) const
    {
        return rootParams_[eIdxGlobal].surface;
    }

    Scalar Kr(unsigned int eIdxGlobal) const
    {
        return rootParams_[eIdxGlobal].radialPerm;
    }

    Scalar rootOrder(unsigned int eIdxGlobal) const
    {
        return rootParams_[eIdxGlobal].order;
    }

    Scalar rootBranch(unsigned int eIdxGlobal) const
    {
        return rootParams_[eIdxGlobal].branchId;
    }

    Scalar rootMass(unsigned int eIdxGlobal) const
    {
        return rootParams_[eIdxGlobal].mass;
    }

    RootParams& rootParams(const Element &element)
    {
        auto eIdx = gridView_.indexSet().index(element);
        return rootParams_[eIdx];
    }
    const RootParams& rootParams(const Element &element) const
    {
        auto eIdx = gridView_.indexSet().index(element);
        return rootParams_[eIdx];
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the porosity
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    { return 0.4; }

private:
    std::vector<RootParams> rootParams_;
    Scalar Kx_, Kr_;
    GridView gridView_;
};

} // end namespace Dumux

#endif
