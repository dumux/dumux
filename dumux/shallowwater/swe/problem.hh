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
 * \ingroup SweModel
 * \copydoc Dumux::SweProblem
 */
#ifndef DUMUX_SWE_PROBLEM_HH
#define DUMUX_SWE_PROBLEM_HH

#include <dumux/common/fvproblem.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include "model.hh"

namespace Dumux
{


/*!
 * \ingroup SweModel
 * \brief Swe problem base class.
 *
 *
 */
template<class TypeTag>
class SweProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);

        /*!
     * \brief Constructor, passing the spatial parameters
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param spatialParams The spatial parameter class
     */
    SweProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                            std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(fvGridGeometry)
    , spatialParams_(spatialParams)
    {}

    /*!
     * \brief Constructor, constructing the spatial parameters
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    SweProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : SweProblem(fvGridGeometry,
                              std::make_shared<SpatialParams>(fvGridGeometry))
    {}


    Scalar getH(const Element& element)
    {
        return this->asImp_().getH();
    }


    Scalar getU(const Element& element)
    {
        return this->asImp_().getU();
    }

    Scalar getV(const Element& element)
    {
        return this->asImp_().getV();
    }

    Scalar getZ(const Element& element)
    {
        return this->asImp_().getZ();
    }

    Scalar getKs(const Element& element)
    {
        return this->asImp_().getKs();
    }

    Scalar getFrictionH(const Element& element)
    {
        return this->asImp_().getFrictionH();
    }

    Scalar getFrictionUstarH(const Element& element)
    {
        return this->asImp_().getFrictionUstarH();
    }

    Scalar getGravity(const Element& element)
    {
        return this->asImp_().getGravity();
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return *spatialParams_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

    // \}


protected:
    GlobalPosition gravity_;

    // material properties
    std::shared_ptr<SpatialParams> spatialParams_;
};

}

#endif
