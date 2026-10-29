// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters in porous medium problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FV_SPATIAL_PARAMS_HH
#define DUMUX_POROUS_MEDIUM_FV_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
template<class GlobalPosition, class SolidSystem>
struct hasInertVolumeFractionAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.template inertVolumeFractionAtPos<SolidSystem>(std::declval<GlobalPosition>(), 0))
    {}
};

template<class GlobalPosition>
struct hasPorosityAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.porosityAtPos(std::declval<GlobalPosition>()))
    {}
};
} // end namespace Detail
#endif

/*!
 * \ingroup Core
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of porous-medium problems.
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVPorousMediumSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    FVPorousMediumSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

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

        return this->asImp_().porosityAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the porosity at an interpolation point.
     * \note There is no default. A local dof on an element face belongs to the elements on both sides,
     *       so a porosity given as a function of the position is ambiguous where it jumps.
     *       Overload this function using the position of the interpolation point for a smooth
     *       porosity and the element of the element discretization for a porosity that is
     *       constant per element.
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    template<class ElementDiscretization, class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization& elemDisc,
                    const IpData& ipData,
                    const ElementSolution& elemSol) const
    {
        static_assert(AlwaysFalse<IpData>::value, " \n\n"
        "   Your spatial params class has to overload\n\n"
        "         template<class ElementDiscretization, class IpData, class ElementSolution>\n"
        "         Scalar porosity(const ElementDiscretization& elemDisc,\n"
        "                         const IpData& ipData,\n"
        "                         const ElementSolution& elemSol) const\n\n"
        "   returning e.g. porosityAtPos(ipData.global()) for a smooth porosity\n"
        "   or a value determined by elemDisc.element() for a porosity that is constant per element.\n\n");
        return 0.0;
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
     * \note this overload is enable if there is only one inert solid component and the
     *       user didn't choose to implement a inertVolumeFractionAtPos overload.
     *       It then forwards to the simpler porosity interface.
     *       With more than one solid components or active solid components (i.e. dissolution)
     *       please overload the more general inertVolumeFraction/inertVolumeFractionAtPos interface.
     */
    template<class SolidSystem, class ElementSolution,
             typename std::enable_if_t<SolidSystem::isInert()
                                       && SolidSystem::numInertComponents == 1
                                       && !decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())(std::declval<Implementation>()))::value,
                                       int> = 0>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        return 1.0 - this->asImp_().porosity(element, scv, elemSol);
    }

    //! \copydoc inertVolumeFraction
    template<class SolidSystem, class IpData, class ElementSolution>
        requires (SolidSystem::isInert()
                  && SolidSystem::numInertComponents == 1
                  && !decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())(std::declval<Implementation>()))::value)
    Scalar inertVolumeFraction(const FVElementGeometry& fvGeometry,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        return 1.0 - this->asImp_().porosity(fvGeometry, ipData, elemSol);
    }

    // specialization if there are no inert components at all
    template<class SolidSystem, class ElementSolution,
             typename std::enable_if_t<SolidSystem::numInertComponents == 0, int> = 0>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        return 0.0;
    }

    //! \copydoc inertVolumeFraction
    template<class SolidSystem, class IpData, class ElementSolution>
        requires (SolidSystem::numInertComponents == 0)
    Scalar inertVolumeFraction(const FVElementGeometry& fvGeometry,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        return 0.0;
    }

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
                                        ),
                                        int> = 0>
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

        return this->asImp_().template inertVolumeFractionAtPos<SolidSystem>(scv.center(), compIdx);
    }

    //! \copydoc inertVolumeFraction
    template<class SolidSystem, class IpData, class ElementSolution>
        requires ((SolidSystem::numInertComponents > 1)
                  || ((SolidSystem::numInertComponents > 0)
                      && (!SolidSystem::isInert()
                          || decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())
                                  (std::declval<Implementation>()))::value)))
    Scalar inertVolumeFraction(const FVElementGeometry& fvGeometry,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    { return this->asImp_().template inertVolumeFractionAtPos<SolidSystem>(ipData.global(), compIdx); }
};

} // namespace Dumux

#endif
