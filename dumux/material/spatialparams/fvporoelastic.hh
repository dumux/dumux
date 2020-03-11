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
 * \brief The base class for spatial parameters of poro-elastic geomechanical problems
 */
#ifndef DUMUX_GEOMECHANICS_POROELASTIC_FV_SPATIAL_PARAMS_HH
#define DUMUX_GEOMECHANICS_POROELASTIC_FV_SPATIAL_PARAMS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class has a lameParamsAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition>
struct hasLameParamsAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.lameParamsAtPos(std::declval<GlobalPosition>()))
    {}
};

// helper struct detecting if the user-defined spatial params class has a reactiveVolumeFractionAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition, class SolidSystem>
struct hasReactiveVolumeFractionAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.template reactiveVolumeFractionAtPos<SolidSystem>(std::declval<GlobalPosition>(), 0))
    {}
};

// helper struct detecting if the user-defined spatial params class has a biotCoefficientAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition>
struct hasBiotCoeffAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.biotCoefficientAtPos(std::declval<GlobalPosition>()))
    {}
};

} // end namespace Detail
#endif

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of poro-elastic geomechanical problems
 */
template<class Scalar, class GridGeometry, class Implementation>
class FVSpatialParamsPoroElastic
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum { dimWorld = GridView::dimensionworld };

public:
    //! The constructor
    FVSpatialParamsPoroElastic(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , gravity_(0.0)
    {
        const bool enableGravity = getParam<bool>("Problem.EnableGravity");
        if (enableGravity)
            gravity_[dimWorld-1]  = -9.81;
    }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * The default behaviour is a constant gravity vector;
     * if the <tt>Problem.EnableGravity</tt> parameter is true,
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     *
     * \param pos the spatial position at which to evaulate the gravity vector
     */
    const GlobalPosition& gravity(const GlobalPosition &pos) const
    { return gravity_; }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     * \note this can only be used for solids with one inert component
     *       (see inertVolumeFraction for the more general interface)
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static_assert(decltype(isValid(Detail::hasPorosityAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         Scalar porosityAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         Scalar porosity(const Element& element,\n"
        "                         const SubControlVolume& scv,\n"
        "                         const ElementSolution& elemSol) const\n\n");

        return asImp_().porosityAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     *
     * \note this overload is enabled if there is only one inert solid component and the
     *       user didn't choose to implement an inertVolumeFractionAtPos overload.
     *       It then forwards to the simpler porosity interface.
     *       With more than one solid components or active solid components (i.e. dissolution)
     *       please overload the more general inertVolumeFraction/inertVolumeFractionAtPos interface.
     */
    template<class SolidSystem, class ElementSolution,
             typename std::enable_if_t<SolidSystem::isInert()
                                       && SolidSystem::numInertComponents == 1
                                       && !decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())
                                                           (std::declval<Implementation>()))::value, int> = 0>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    { return 1.0 - asImp_().porosity(element, scv, elemSol); }

    // specialization if there are no inert components at all
    template<class SolidSystem, class ElementSolution, typename std::enable_if_t<SolidSystem::numInertComponents == 0, int> = 0>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    { return 0.0; }

    // the more general interface forwarding to inertVolumeFractionAtPos
    template<class SolidSystem, class ElementSolution,
             typename std::enable_if_t<(SolidSystem::numInertComponents > 1) ||
                                       (
                                            (SolidSystem::numInertComponents > 0) &&
                                            (
                                                !SolidSystem::isInert()
                                                || decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())
                                                                   (std::declval<Implementation>()))::value
                                            )
                                        ), int> = 0>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        static_assert(decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         template<class SolidSystem>\n"
        "         Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const\n\n"
        "   or overload this function\n\n"
        "         template<class SolidSystem, class ElementSolution>\n"
        "         Scalar inertVolumeFraction(const Element& element,\n"
        "                                    const SubControlVolume& scv,\n"
        "                                    const ElementSolution& elemSol,\n"
        "                                    int compIdx) const\n\n");

        return asImp_().template inertVolumeFractionAtPos<SolidSystem>(scv.center(), compIdx);
    }

    /*!
     * \brief Function for defining the solid volume fraction of a solid
     *        component that takes part in some sort of reaction.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     *
     * \note This overload is enabled if there are only inert solid components
     *       and the user did not choose to implement a reactiveVolumeFractionAtPos
     *       function. The reactive volume fraction is zero in this case.
     */
    template<class SolidSystem, class ElementSolution,
             std::enable_if_t< SolidSystem::isInert() &&
                               !decltype(isValid(Detail::hasReactiveVolumeFractionAtPos<GlobalPosition, SolidSystem>())
                                                (std::declval<Implementation>()))::value, int > = 0 >
    Scalar reactiveVolumeFraction(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol,
                                  int compIdx) const
    { return 0.0; }

    //! overload for the case of reactive solids or user-provided overload
    template<class SolidSystem, class ElementSolution,
             std::enable_if_t< !SolidSystem::isInert() ||
                               decltype(isValid(Detail::hasReactiveVolumeFractionAtPos<GlobalPosition, SolidSystem>())
                                                (std::declval<Implementation>()))::value, int > = 0 >
    Scalar reactiveVolumeFraction(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol,
                                  int compIdx) const
    {
        static_assert(decltype(isValid(Detail::hasReactiveVolumeFractionAtPos<GlobalPosition, SolidSystem>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         template<class SolidSystem>\n"
        "         Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const\n\n"
        "   or overload this function\n\n"
        "         template<class SolidSystem, class ElementSolution>\n"
        "         Scalar inertVolumeFraction(const Element& element,\n"
        "                                    const SubControlVolume& scv,\n"
        "                                    const ElementSolution& elemSol,\n"
        "                                    int compIdx) const\n\n");

        return asImp_().template reactiveVolumeFractionAtPos<SolidSystem>(scv.center(), compIdx);
    }

    /*!
     * \brief Define the Lame parameters
     * \note  These are possibly solution dependent and are evaluated
     *        for an integration point inside the element. Therefore,
     *        a flux variables cache object is passed to this function
     *        containing data on shape functions at the integration point.
     *
     * \param element The current element
     * \param fvGeometry The local finite volume geometry
     * \param elemVolVars Primary/Secondary variables inside the element
     * \param fluxVarsCache Contains data on shape functions at the integration point
     * \return lame parameters
     */
    template<class ElemVolVars, class FluxVarsCache>
    decltype(auto) lameParams(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElemVolVars& elemVolVars,
                              const FluxVarsCache& fluxVarsCache) const
    {
        static_assert(decltype(isValid(Detail::hasLameParamsAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         const LameParams& lameParams(const Element& element,\n"
        "                                      const FVElementGeometry& fvGeometry,\n"
        "                                      const ElemVolVars& elemVolVars,\n"
        "                                      const FluxVarsCache& fluxVarsCache) const\n\n");

        return asImp_().lameParamsAtPos(fluxVarsCache.ipGlobal());
    }

    /*!
     * \brief Returns the Biot coefficient in an element
     * \note  This is possibly solution dependent and is evaluated
     *        for an integration point inside the element. Therefore,
     *        a flux variables cache object is passed to this function
     *        containing data on shape functions at the integration point.
     *
     * \param element The current element
     * \param fvGeometry The local finite volume geometry
     * \param elemVolVars Primary/Secondary variables inside the element
     * \param fluxVarsCache Contains data on shape functions at the integration point
     * \return Biot coefficient
     */
    template<class ElemVolVars, class FluxVarsCache>
    Scalar biotCoefficient(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElemVolVars& elemVolVars,
                           const FluxVarsCache& fluxVarsCache) const
    {
        static_assert(decltype(isValid(Detail::hasBiotCoeffAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const LameParams& biotCoefficientAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         const LameParams& biotCoefficient(const Element& element,\n"
        "                                      const FVElementGeometry& fvGeometry,\n"
        "                                      const ElemVolVars& elemVolVars,\n"
        "                                      const FluxVarsCache& fluxVarsCache) const\n\n");

        return asImp_().biotCoefficientAtPos(fluxVarsCache.ipGlobal());
    }

    //! The finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    GlobalPosition gravity_; //!< The gravity vector
};
} // end namespace Dumux
#endif
