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
 * \ingroup SweTests
 * \brief The spatial parameters class for the test problem using the
 *        shallow water model
 */
#ifndef DUMUX_SWE_TEST_SPATIALPARAMS_HH
#define DUMUX_SWE_TEST_SPATIALPARAMS_HH


namespace Dumux
{

//forward declaration
template<class TypeTag>
class SweTestSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(SweTestSpatialParams);
}

/*!
 * \ingroup SweTests
 * \brief The spatial parameters class for the test problem using the
 *        shallow water model
 */
template<class TypeTag>
class SweTestSpatialParams //: public FVSpatialParamsSwe<TypeTag>
{
    //using ParentType = FVSpatialParamsSwe<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexSet = typename GridView::IndexSet;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    SweTestSpatialParams(const Problem& problem)
    {}

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*! \brief Define the friction parameter ks.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    * \param elemSol The solution at the dofs connected to the element.
    * \return the material parameters object
    */
    Scalar ks(const Element& element,
             const SubControlVolume& scv,
             const ElementSolutionVector& elemSol) const
    {
        // Get the global index of the element
        //const auto eIdx = this->problem().fvGridGeometry().elementMapper().index(element);
        return 0.0;
    }

    /*! \brief Define the gravitation.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    * \param elemSol The solution at the dofs connected to the element.
    * \return the material parameters object
    */
    Scalar grav() const
    {
        return grav_;
    }

    /*! \brief Define the friciton law.
    *
    *   0 = no friction
    *   1 = Manning
    *   2 = Chezy
    *   3 = Nikuradse
    *
    *
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    * \param elemSol The solution at the dofs connected to the element.
    * \return the material parameters object
    */
    int frictionlaw() const
    {
        return 0;
    }


    Scalar bottom() const
    {
        return 0.0;
    }


private:

//    const IndexSet& indexSet_;
    static constexpr Scalar eps_ = 1.5e-7;
    static constexpr Scalar grav_ = 9.81;
};

} // end namespace

#endif
