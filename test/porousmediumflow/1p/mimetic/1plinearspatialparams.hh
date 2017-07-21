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
#ifndef DUMUX_1P_LINEAR_SPATIALPARAMS_HH
#define DUMUX_1P_LINEAR_SPATIALPARAMS_HH

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
class OnePLinearSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexSet = typename GridView::IndexSet;
    using ScalarVector = std::vector<Scalar>;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    OnePLinearSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView),permOne_(0),permTwo_(0)
    {
            permOne_[0][0] = 3.0;
            permOne_[1][1] = 3.0;
            permOne_[0][1] = permOne_[1][0] = 1.0;

            if(dimWorld == 3)
                permOne_[2][2] = 1.0;

            permTwo_[0][0] = 10.0;
            permTwo_[1][1] = 10.0;
            permTwo_[0][1] = permTwo_[1][0] = 3.0;

            if(dimWorld == 3)
                permTwo_[2][2] = 1.0;

            testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, TestCase);
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
        const GlobalPosition &globalPos = scv.center();

        PermeabilityType perm(0);

        if(testCase_ == 1)
        {
            if(globalPos[0]<0.6+1.0e-8)
                perm = permOne_;
            else
                perm = permTwo_;
        }
        else
        {
            perm[0][0] = perm[1][1] = 1.0;
            perm[dimWorld-1][dimWorld-1] = 1.0e-3;
        }

        return perm;
    }


    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    int testCase_;
    DimWorldMatrix permOne_;
    DimWorldMatrix permTwo_;
};
} // end namespace
#endif
