// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of problems using the
 *        fv method.
 */
#ifndef DUMUX_SEQUENTIAL_FV_SPATIAL_PARAMS_HH
#define DUMUX_SEQUENTIAL_FV_SPATIAL_PARAMS_HH

#include <dumux/common/properties.hh>
#include "sequentialfv1p.hh"

namespace Dumux {
namespace Properties
{ template<class TypeTag, class MyTypeTag>
struct  MaterialLaw  { using type = UndefinedProperty; }; }

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of a multi-phase problem using the
 *        fv method.
 */
template<class TypeTag>
class SequentialFVSpatialParams: public SequentialFVSpatialParamsOneP<TypeTag>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Implementation = GetPropType<TypeTag, Properties::SpatialParams>;

    enum
    {
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    SequentialFVSpatialParams(const Problem& problem)
    :SequentialFVSpatialParamsOneP<TypeTag>(problem)
    {
    }


    // make sure we get a deprecation warning (remove this after release 3.3)
    template<class ElementSolution, class SubControlVolume>
    [[deprecated("Use the new style material laws. Old material laws and this interface will no longer be supported after release 3.3")]]
    decltype(auto) materialLawParamsDeprecated(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        return this->asImp_().materialLawParams(element);
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \return the material parameters object
     * \param element The element
     */
    const auto& materialLawParams(const Element &element) const
    {
        return asImp_().materialLawParamsAtPos(element.geometry().center());
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \return the material parameters object
     * \param globalPos The position of the center of the element
     */
    const auto& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a materialLawParamsAtPos() method.");
    }

private:
    Implementation &asImp_()
    {
        return *static_cast<Implementation*> (this);
    }

    const Implementation &asImp_() const
    {
        return *static_cast<const Implementation*> (this);
    }
};

} // namespace Dumux

#endif
