// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters in multi-phase porous-medium-flow problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_MP_HH
#define DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_MP_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include "fvspatialparams.hh"

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class
// has a fluidMatrixInteractionAtPos function
template<class GlobalPosition>
struct hasFluidMatrixInteractionAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.fluidMatrixInteractionAtPos(std::declval<GlobalPosition>()))
    {}
};
} // end namespace Detail
#endif

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of multi-phase problems
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVPorousMediumFlowSpatialParamsMP
: public FVPorousMediumFlowSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    FVPorousMediumFlowSpatialParamsMP(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class ElementSolution>
    decltype(auto) fluidMatrixInteraction(const Element& element,
                                          const SubControlVolume& scv,
                                          const ElementSolution& elemSol) const
    {
        static_assert(decltype(isValid(Detail::hasFluidMatrixInteractionAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         auto fluidMatrixInteraction(const Element& element,\n"
        "                                     const SubControlVolume& scv,\n"
        "                                     const ElementSolution& elemSol) const\n\n");

        return this->asImp_().fluidMatrixInteractionAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.)
     *        at an interpolation point.
     * \note There is no default. A local dof on an element face belongs to the elements on both sides,
     *       so parameters given as a function of the position are ambiguous where they jump.
     *       Overload this function using the position of the interpolation point for smooth
     *       parameters and the element of the element discretization for parameters that are
     *       constant per element.
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class ElementDiscretization, class IpData, class ElementSolution>
    void fluidMatrixInteraction(const ElementDiscretization& elemDisc,
                                const IpData& ipData,
                                const ElementSolution& elemSol) const
    {
        static_assert(AlwaysFalse<IpData>::value, " \n\n"
        "   Your spatial params class has to overload\n\n"
        "         template<class ElementDiscretization, class IpData, class ElementSolution>\n"
        "         auto fluidMatrixInteraction(const ElementDiscretization& elemDisc,\n"
        "                                     const IpData& ipData,\n"
        "                                     const ElementSolution& elemSol) const\n\n"
        "   returning e.g. fluidMatrixInteractionAtPos(ipData.global()) for smooth parameters\n"
        "   or a value determined by elemDisc.element() for parameters that are constant per element.\n\n");
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
     * \brief Function for defining which phase is to be considered as the wetting phase
     *        at an interpolation point.
     * \note There is no default, for the same reason as for the parameters of the constitutive relationships
     *       at an interpolation point.
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class FluidSystem, class ElementDiscretization, class IpData, class ElementSolution>
    int wettingPhase(const ElementDiscretization& elemDisc,
                     const IpData& ipData,
                     const ElementSolution& elemSol) const
    {
        static_assert(AlwaysFalse<IpData>::value, " \n\n"
        "   Your spatial params class has to overload\n\n"
        "         template<class FluidSystem, class ElementDiscretization, class IpData, class ElementSolution>\n"
        "         int wettingPhase(const ElementDiscretization& elemDisc,\n"
        "                          const IpData& ipData,\n"
        "                          const ElementSolution& elemSol) const\n\n"
        "   returning e.g. wettingPhaseAtPos<FluidSystem>(ipData.global()) if the wetting phase is the same\n"
        "   everywhere or a value determined by elemDisc.element() if it is constant per element.\n\n");
        return 0;
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
