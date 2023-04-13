// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters in porous-medium-flow problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_HH
#define DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/fvporousmediumspatialparams.hh>
#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined spatial params class has a permeabilityAtPos function
template<class GlobalPosition>
struct hasPermeabilityAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.permeabilityAtPos(std::declval<GlobalPosition>()))
    {}
};
} // end namespace Detail
#endif

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of porous-medium-flow problems.
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVPorousMediumFlowSpatialParams
: public FVPorousMediumSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumSpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    FVPorousMediumFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
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
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    template<class ElementSolution>
    decltype(auto) permeability(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        static_assert(decltype(isValid(Detail::hasPermeabilityAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const PermeabilityType& permeabilityAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         const PermeabilityType& permeability(const Element& element,\n"
        "                                              const SubControlVolume& scv,\n"
        "                                              const ElementSolution& elemSol) const\n\n");

        return this->asImp_().permeabilityAtPos(scv.center());
    }

    /*!
     * \brief If the permeability should be evaluated directly at the scvf integration point
     *        (for convergence tests with analytical and continuous perm functions) or is evaluated
     *        at the scvs (for permeability fields with discontinuities) -> default
     */
    static constexpr bool evaluatePermeabilityAtScvfIP()
    { return false; }

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
     * \param scvf The sub-control volume face where the
     *           intrinsic velocity ought to be calculated.
     */
    Scalar forchCoeff(const SubControlVolumeFace &scvf) const
    {
        return forchCoeffDefault_;
    }

private:
    Scalar forchCoeffDefault_;
};

} // namespace Dumux

#endif
