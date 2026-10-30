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
#ifndef DUMUX_COMMON_POROUS_MEDIUM_SPATIAL_PARAMS_HH
#define DUMUX_COMMON_POROUS_MEDIUM_SPATIAL_PARAMS_HH

#include <utility>

#include <dumux/common/spatialparams.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux::Experimental {

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
} // end namespace Detail
#endif

/*!
 * \ingroup Core
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of porous-medium problems.
 */
template<class GridDiscretization, class Scalar, class Implementation>
class PorousMediumSpatialParams
: public SpatialParams<GridDiscretization, Scalar, Implementation>
{
    using ParentType = SpatialParams<GridDiscretization, Scalar, Implementation>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GlobalPosition = typename GridDiscretization::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    using ParentType::ParentType;

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
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
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
    template<class SolidSystem, class IpData, class ElementSolution>
        requires (SolidSystem::isInert()
                  && SolidSystem::numInertComponents == 1
                  && !decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())(std::declval<Implementation>()))::value)
    Scalar inertVolumeFraction(const ElementDiscretization& elemDisc,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        return 1.0 - this->asImp_().porosity(elemDisc, ipData, elemSol);
    }

    //! \copydoc inertVolumeFraction
    template<class SolidSystem, class IpData, class ElementSolution>
        requires (SolidSystem::numInertComponents == 0)
    Scalar inertVolumeFraction(const ElementDiscretization& elemDisc,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        return 0.0;
    }

    //! \copydoc inertVolumeFraction
    template<class SolidSystem, class IpData, class ElementSolution>
        requires ((SolidSystem::numInertComponents > 1)
                  || ((SolidSystem::numInertComponents > 0)
                      && (!SolidSystem::isInert()
                          || decltype(isValid(Detail::hasInertVolumeFractionAtPos<GlobalPosition, SolidSystem>())
                                  (std::declval<Implementation>()))::value)))
    Scalar inertVolumeFraction(const ElementDiscretization& elemDisc,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    { return this->asImp_().template inertVolumeFractionAtPos<SolidSystem>(ipData.global(), compIdx); }
};

} // end namespace Dumux::Experimental

#endif
