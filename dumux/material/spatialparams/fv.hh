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
 * \brief The base class for spatial parameters of multi-phase problems
 * using a fully implicit discretization method.
 */
#ifndef DUMUX_FV_SPATIAL_PARAMS_HH
#define DUMUX_FV_SPATIAL_PARAMS_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include "fv1p.hh"

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class has a materialLawParamsAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition>
struct hasMaterialLawParamsAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.materialLawParamsAtPos(std::declval<GlobalPosition>()))
    {}
};
} // end namespace Detail
#endif


/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of multi-phase problems
 * using a fully implicit discretization method.
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVSpatialParams : public FVSpatialParamsOneP<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    FVSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    // make sure we get a deprecation warning (remove this after release 3.3)
    template<class ElementSolution>
    [[deprecated("Use the new style material laws. Old material laws and this interface will no longer be supported after release 3.3")]]
    decltype(auto) materialLawParamsDeprecated(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        return this->asImp_().materialLawParams(element, scv, elemSol);
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the material parameters object
     */
    template<class ElementSolution>
    [[deprecated("Use the new style material laws. Old material laws and this interface will no longer be supported after release 3.3")]]
    decltype(auto) materialLawParams(const Element& element,
                                     const SubControlVolume& scv,
                                     const ElementSolution& elemSol) const
    {
        static_assert(decltype(isValid(Detail::hasMaterialLawParamsAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         const MaterialLawParams& materialLawParams(const Element& element,\n"
        "                                                    const SubControlVolume& scv,\n"
        "                                                    const ElementSolution& elemSol) const\n\n");

        return this->asImp_().materialLawParamsAtPos(scv.center());
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the wetting phase index
     */
    template<class FluidSystem, class ElementSolution>
    int wettingPhase(const Element& element,
                     const SubControlVolume& scv,
                     const ElementSolution& elemSol) const
    {
        return this->asImp_().template wettingPhaseAtPos<FluidSystem>(scv.center());
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a wettingPhaseAtPos() method.");
    }
};

} // namespace Dumux

#endif
