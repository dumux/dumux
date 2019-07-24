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
 * \brief The base class for spatial parameters of linear elastic geomechanical problems.
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_GEOMECHANICS_ELASTIC_SPATIAL_PARAMS_HH

#include <memory>

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/method.hh>

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

} // end namespace Detail
#endif

//! Forward declaration of the Implementation
template<class Scalar, class GridGeometry, class Implementation, DiscretizationMethod dm>
class SpatialParamsElasticImpl;

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of linear elastic geomechanical problems.
 *        Derive from this class in your spatial parameters implementation and provide overloads
 *        of the respective functions for defining the lame parameters and the inert volume fraction,
 *        where the latter is only required if there are more than one inert components
 * \note Currently, the elastic model works with the box and finite element schemes. The interfaces
 *       to define the parameters are different for the two methods (unless you restrict yourself to
 *       the only position-dependent ...atPos() variant). In your implementation you must overload the
 *       version that fits your choice of discretization method.
 */
template<class Scalar, class GridGeometry, class Implementation>
using SpatialParamsElastic = SpatialParamsElasticImpl<Scalar, GridGeometry, Implementation, GridGeometry::discMethod>;

/*!
 * \ingroup SpatialParameters
 * \brief The common base class for the discretization method-specific implementations.
 */
template<class Scalar, class GridGeometry, class Implementation>
class SpatialParamsElasticBase
{
    using ElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    //! The constructor
    SpatialParamsElasticBase(std::shared_ptr<const GridGeometry> gridGeometry)
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
    const GlobalPosition& gravity(const GlobalPosition& pos) const
    { return gravity_; }

    /*!
     * \brief Return reference to the grid geometry.
     */
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

protected:
    //! return reference to the implementation
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    //! return const reference to the implementation
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    GlobalPosition gravity_; //!< The gravity vector
};

/*!
 * \ingroup SpatialParameters
 * \brief Specialization of the implementation for the box scheme.
 */
template<class Scalar, class GridGeometry, class Implementation>
class SpatialParamsElasticImpl<Scalar, GridGeometry, Implementation, DiscretizationMethod::box>
: public SpatialParamsElasticBase<Scalar, GridGeometry, Implementation>
{
    using ParentType = SpatialParamsElasticBase<Scalar, GridGeometry, Implementation>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using FVElementGeometry = typename GridGeometry::LocalView;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Pull up the parents' constructor
    using ParentType::ParentType;

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     */
    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        static_assert(SolidSystem::isInert(), "Elastic model can only be used with inert solid systems");

        // when there is only one component, the volume fraction must be one
        if (SolidSystem::numInertComponents == 1)
            return 1.0;

        // otherwise we require the user to define the solid composition
        return this->asImp_().template inertVolumeFractionAtPos<SolidSystem>(scv.center(), compIdx);
    }

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param globalPos The global position
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     */
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "The spatial parameters do not provide inertVolumeFractionAtPos() method."); }

    /*!
     * \brief Define the Lame parameters
     * \note  These are possibly solution dependent and are evaluated
     *        for an integration point inside the element. Therefore,
     *        a flux variables cache object is passed to this function
     *        containing data on shape functions at the integration point.
     *
     * \param element The current element
     * \param fvGeometry The local finite volume geometry
     * \param elemVolVars element volume variables of the element (primary/secondary variables)
     * \param fluxVarsCache Contains data on shape functions at the integration point
     * \return lame parameters
     */
    template<class ElementVolumeVariables, class FluxVarsCache>
    decltype(auto) lameParams(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
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

        return this->asImp_().lameParamsAtPos(fluxVarsCache.ipGlobal());
    }

    //! Overload for backwards compatibility
    [[deprecated("Use gridGeometry() instead.")]]
    const GridGeometry& fvGridGeometry() const
    { return this->gridGeometry(); }
};

/*!
 * \ingroup SpatialParameters
 * \brief Specialization of the implementation for finite element schemes.
 */
template<class Scalar, class GridGeometry, class Implementation>
class SpatialParamsElasticImpl<Scalar, GridGeometry, Implementation, DiscretizationMethod::fem>
: public SpatialParamsElasticBase<Scalar, GridGeometry, Implementation>
{
    using ParentType = SpatialParamsElasticBase<Scalar, GridGeometry, Implementation>;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Pull up the parents' constructor
    using ParentType::ParentType;

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param ipData The shape function values/gradients evaluated at an integration point.
     * \param elemSol The solution at the dofs connected to the element.
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     */
    template<class SolidSystem, class IpData, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const IpData& ipData,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        static_assert(SolidSystem::isInert(), "Elastic model can only be used with inert solid systems");

        // when there is only one component, the volume fraction must be one
        if (SolidSystem::numInertComponents == 1)
            return 1.0;

        // otherwise we require the user to define the solid composition
        return this->asImp_().template inertVolumeFractionAtPos<SolidSystem>(ipData.ipGlobal(), compIdx);
    }

    /*!
     * \brief Function for defining the solid volume fraction.
     *        That is possibly solution dependent.
     *
     * \param globalPos The global position
     * \param compIdx The solid component index
     * \return the volume fraction of the inert solid component with index compIdx
     */
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "The spatial parameters do not provide inertVolumeFractionAtPos() method."); }

    /*!
     * \brief Define the Lame parameters
     * \note  These are possibly solution dependent and are evaluated
     *        for an integration point inside the element. Therefore,
     *        an ipData object is passed to this function containing
     *        data on shape functions at the integration point.
     *
     * \param element The current element
     * \param feGeometry The local finite element geometry
     * \param elemSol The element solution vector
     * \param ipData Contains data on shape functions at the integration point
     * \param secVars Contains the primary/secondary variables at integration point
     * \return lame parameters
     */
    template<class ElementSolution, class IpData, class SecondaryVariables>
    decltype(auto) lameParams(const Element& element,
                              const FEElementGeometry& feGeometry,
                              const ElementSolution& elemSol,
                              const IpData& ipData,
                              const SecondaryVariables& secVars) const
    {
        static_assert(decltype(isValid(Detail::hasLameParamsAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         const LameParams& lameParams(const Element& element,\n"
        "                                      const FVElementGeometry& fvGeometry,\n"
        "                                      const ElemIpDataVolVars& ipData,\n"
        "                                      const SecondaryVariables& secVars) const\n\n");

        return this->asImp_().lameParamsAtPos(ipData.ipGlobal());
    }
};
} // end namespace Dumuxs
#endif
