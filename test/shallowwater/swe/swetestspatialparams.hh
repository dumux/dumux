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
    { return 1.0; }

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
    * \return the friction law of the object
    */
    int frictionlaw() const
    {
        return 0;
    }

    /*! \brief Define the bottom (for SWEs without sediment transport
    *          this should be z.
    *
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    * \return the material parameters object
    */
    Scalar bottom(const Element& element,
             const SubControlVolume& scv) const
    {
        auto someInitSol = element.geometry().center();
        /*if ((11.0 < someInitSol[0])&&(someInitSol[0] < 14.0)){
            return 3.1;
        }*/
        return 0.0;
    }

    /*! \brief Set all element data for SWEs and set them.
    *
    *
    * \param Gridview gv
    * \param elementdata map
    * \return the material parameters object
    */
    void setElementdata(auto elementdata)
    {
        //set all data
        if ( elementdata.find("z") != elementdata.end() ){
            eleZ_ = elementdata["z"];
        }else{
            //TODO thorw an error
            std::cout << "Error elementdata missing z" << std::endl;
        }

        //optional variables
        auto nelem = eleZ_.size();

        if (elementdata.find("ks") != elementdata.end() ){
            eleKs_ = elementdata["ks"];
        }else{
            eleKs_.assign(nelem, 0.0);
        }

        if (elementdata.find("vegDp") != elementdata.end() ){
            eleVegDp_ = elementdata["vegDp"];
        }else{
            eleVegDp_.assign(nelem, 0.0);
        }

        if (elementdata.find("vegAp") != elementdata.end() ){
            eleVegAp_ = elementdata["vegAp"];
        }else{
            eleVegAp_.assign(nelem, 0.0);
        }

        if (elementdata.find("vegHp") != elementdata.end() ){
            eleVegHp_ = elementdata["vegHp"];
        }else{
            eleVegHp_.assign(nelem, 0.0);
        }

        if (elementdata.find("zoneId") != elementdata.end() ){
            eleZoneId_ = elementdata["zoneId"];
        }else{
            eleZoneId_.assign(nelem, 0.0);
        }

        if (elementdata.find("defaultId") != elementdata.end() ){
            eleDefaultId_ = elementdata["defaultId"];
        }else{
            eleDefaultId_.assign(nelem, 0.0);
        }
    }


private:

   //const IndexSet& indexSet_;
    static constexpr Scalar eps_ = 1.5e-7;
    static constexpr Scalar grav_ = 9.81;

    // element data
    std::vector<double> eleZ_;
    std::vector<double> eleKs_;
    std::vector<double> eleVegDp_;
    std::vector<double> eleVegAp_;
    std::vector<double> eleVegHp_;
    std::vector<double> eleCflLentgth_;
    std::vector<double> eleDefaultId_;
    std::vector<double> eleZoneId_; //Do we need this?
};

} // end namespace

#endif
