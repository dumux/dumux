// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters in porous-medium-flow problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_SPATIAL_PARAMS_HH
#define DUMUX_POROUS_MEDIUM_FLOW_SPATIAL_PARAMS_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/porousmediumspatialparams.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of porous-medium-flow problems.
 */
template<class GridDiscretization, class Scalar, class Implementation>
class PorousMediumFlowSpatialParams
: public PorousMediumSpatialParams<GridDiscretization, Scalar, Implementation>
{
    using ParentType = PorousMediumSpatialParams<GridDiscretization, Scalar, Implementation>;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    PorousMediumFlowSpatialParams(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    {
        /* \brief default forchheimer coefficient
         * Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90 \cite ward1964 .
         *        Actually the Forchheimer coefficient is also a function of the dimensions of the
         *        porous medium. Taking it as a constant is only a first approximation
         *        (Nield, Bejan, Convection in porous media, 2006, p. 10 \cite nield2006 )
         */
        forchCoeffDefault_ = getParam<Scalar>("SpatialParams.ForchCoeff", 0.55);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$ at an
     *        interpolation point.
     * \note There is no default. A local dof on an element face belongs to the elements on both sides,
     *       so a permeability given as a function of the position is ambiguous where it jumps.
     *       Overload this function using the position of the interpolation point for a smooth
     *       permeability and the element of the element discretization for a permeability that is
     *       constant per element.
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     * \param elemSol The solution at the dofs connected to the element.
     * \return the intrinsic permeability
     */
    template<class ElementDiscretization, class IpData, class ElementSolution>
    void permeability(const ElementDiscretization& elemDisc,
                      const IpData& ipData,
                      const ElementSolution& elemSol) const
    {
        static_assert(AlwaysFalse<IpData>::value, " \n\n"
        "   Your spatial params class has to overload\n\n"
        "         template<class ElementDiscretization, class IpData, class ElementSolution>\n"
        "         PermeabilityType permeability(const ElementDiscretization& elemDisc,\n"
        "                                       const IpData& ipData,\n"
        "                                       const ElementSolution& elemSol) const\n\n"
        "   returning e.g. permeabilityAtPos(ipData.global()) for a smooth permeability\n"
        "   or a value determined by elemDisc.element() for a permeability that is constant per element.\n\n");
    }

    /*!
     * \brief Function for defining the Beavers-Joseph coefficient for multidomain
     *        problems\f$\mathrm{[-]}\f$.
     *
     * \return Beavers-Joseph coefficient \f$\mathrm{[-]}\f$
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide a beaversJosephCoeffAtPos() method.");
    }

    /*!
     * \brief Apply the Forchheimer coefficient for inertial forces
     *        calculation.
     * \param ipData The interpolation point data where the
     *           intrinsic velocity ought to be calculated.
     */
    template<class IpData>
    Scalar forchCoeff(const IpData& ipData) const
    { return forchCoeffDefault_; }

private:
    Scalar forchCoeffDefault_;
};

} // end namespace Dumux::Experimental

#endif
