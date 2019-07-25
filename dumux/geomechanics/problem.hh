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
 * \ingroup Geomechanics
 * \brief Base class for all geomechanical problems
 */
#ifndef DUMUX_GEOMECHANICS_PROBLEM_HH
#define DUMUX_GEOMECHANICS_PROBLEM_HH

#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/feproblem.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper struct detecting if the user-defined problem class has an effectiveFluidDensityAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition>
struct hasEffFluidDensityAtPos
{
    template<class Problem>
    auto operator()(const Problem& a)
    -> decltype(a.effectiveFluidDensityAtPos(std::declval<GlobalPosition>()))
    {}
};

// helper struct detecting if the user-defined problem class has an effectivePorePressureAtPos function
// for g++ > 5.3, this can be replaced by a lambda
template<class GlobalPosition>
struct hasEffPorePressureAtPos
{
    template<class Problem>
    auto operator()(const Problem& a)
    -> decltype(a.effectivePorePressureAtPos(std::declval<GlobalPosition>()))
    {}
};

} // end namespace Detail
#endif

//! Forward declaration of the implementation
template<class TypeTag, DiscretizationMethod dm>
class GeomechanicsProblemImpl;

/*!
 * \ingroup Geomechanics
 * \brief Base class for all geomechanical problems
 * \note The geomechanics module currently works with the box scheme
 *       and finite element schemes. This alias selects the right
 *       implementation of the problem base class.
 */
template<class TypeTag>
using GeomechanicsProblem
= GeomechanicsProblemImpl<TypeTag, GetPropType<TypeTag, Properties::GridGeometry>::discMethod>;

/*!
 * \ingroup Geomechanics
 * \brief Implementation of the base class for the box scheme
 * \note We require only little additional functionality to the
 *       porous medium flow problem, which is why we inherit from that here.
 */
template<class TypeTag>
class GeomechanicsProblemImpl<TypeTag, DiscretizationMethod::box>
: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int numFP = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();

public:
    //! pull up the constructor of the parent class
    using ParentType::ParentType;

    /*!
     * \brief Returns the effective fluid density within an scv.
     * \note This is only enabled if the model considers fluid phases.
     *
     * \param element The current element
     * \param scv The sub-control volume
     */
    template< int n = numFP, std::enable_if_t<(n > 0), int> = 0 >
    Scalar effectiveFluidDensity(const Element& element,
                                 const SubControlVolume& scv) const
    {
        static_assert(decltype(isValid(Detail::hasEffFluidDensityAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your problem class has to either implement\n\n"
        "         Scalar effectiveFluidDensityAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         Scalar effectiveFluidDensity(const Element& element,\n\
                                               const SubControlVolume& scv) const\n\n");

        return this->asImp_().effectiveFluidDensityAtPos(scv.center());
    }

    /*!
     * \brief Returns the effective pore pressure
     * \note This is only enabled if the model considers fluid phases.
     *       This is possibly solution dependent and is evaluated
     *       for an integration point inside the element. Therefore,
     *       a flux variables cache object is passed to this function
     *       containing data on shape functions at the integration point.
     *
     * \param element The current element
     * \param fvGeometry The local finite volume geometry
     * \param elemVolVars Primary/Secondary variables inside the element
     * \param fluxVarsCache Contains data on shape functions at the integration point
     */
    template< class ElemVolVars, class FluxVarsCache, int n = numFP, std::enable_if_t<(n > 0), int> = 0 >
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElemVolVars& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    {
        static_assert(decltype(isValid(Detail::hasEffPorePressureAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your problem class has to either implement\n\n"
        "         Scalar effectivePorePressureAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         Scalar effectivePorePressure(const Element& element,\n"
        "                                      const FVElementGeometry& fvGeometry,\n"
        "                                      const ElemVolVars& elemVolVars,\n"
        "                                      const FluxVarsCache& fluxVarsCache) const\n\n");

        return this->asImp_().effectivePorePressureAtPos(element.geometry().center());
    }
};

/*!
 * \ingroup Geomechanics
 * \brief Implementation of the base class for finite element schemes
 */
template<class TypeTag>
class GeomechanicsProblemImpl<TypeTag, DiscretizationMethod::fem>
: public FEProblem<TypeTag>
{
    using ParentType = FEProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SecondaryVariables = typename GridVariables::SecondaryVariables;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int numFP = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();

public:
    //! Export spatial parameter type
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    /*!
     * \brief Constructor, passing the spatial parameters.
     *
     * \param gridGeometry The grid geometry
     * \param spatialParams The spatial parameter class
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    GeomechanicsProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<SpatialParams> spatialParams,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , spatialParams_(spatialParams)
    {}

    /*!
     * \brief Constructor, constructing the spatial parameters.
     *
     * \param gridGeometry The grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    GeomechanicsProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : GeomechanicsProblemImpl(gridGeometry,
                              std::make_shared<SpatialParams>(gridGeometry),
                              paramGroup)
    {}

    /*!
     * \name Physical parameters for geomechanical problems
     */
    // \{

    /*!
     * \brief Returns the effective fluid density within an scv.
     * \note This is only enabled if the model considers fluid phases.
     *
     * \param element The current element
     * \param ipData Contains data on shape values/gradients at an integration point
     * \param elemSol The element solution vector
     */
    template< class IpData, class ElementSolution, int n = numFP, std::enable_if_t<(n > 0), int> = 0 >
    Scalar effectiveFluidDensity(const Element& element,
                                 const IpData& ipData,
                                 const ElementSolution& elemSol) const
    {
        static_assert(decltype(isValid(Detail::hasEffFluidDensityAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your problem class has to either implement\n\n"
        "         Scalar effectiveFluidDensityAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         Scalar effectiveFluidDensity(const Element& element,\n\
                                               const SubControlVolume& scv) const\n\n");

        return this->asImp_().effectiveFluidDensityAtPos(ipData.ipGlobal());
    }

    /*!
     * \brief Returns the effective pore pressure
     * \note This is only enabled if the model considers fluid phases.
     *       This is possibly solution dependent and is evaluated
     *       for an integration point inside the element. Therefore,
     *       a flux variables cache object is passed to this function
     *       containing data on shape functions at the integration point.
     *
     * \param element The current element
     * \param feGeometry The finite element geometry
     * \param elemSol The element solution vector
     * \param ipData Shape values/gradients at an integration point
     * \param secVars Primary/secondary variables evaluated at the integration point
     */
    template< class ElementSolution, class IpData, int n = numFP, std::enable_if_t<(n > 0), int> = 0 >
    Scalar effectivePorePressure(const Element& element,
                                 const ElementGeometry& feGeometry,
                                 const ElementSolution& elemSol,
                                 const IpData& ipData,
                                 const SecondaryVariables& secVars) const
    {
        static_assert(decltype(isValid(Detail::hasEffPorePressureAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your problem class has to either implement\n\n"
        "         Scalar effectivePorePressureAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         template<class ElementSolution>\n"
        "         Scalar effectivePorePressure(const Element& element,\n"
        "                                      const FVElementGeometry& fvGeometry,\n"
        "                                      const ElemVolVars& elemVolVars,\n"
        "                                      const FluxVarsCache& fluxVarsCache) const\n\n");

        return this->asImp_().effectivePorePressureAtPos(ipData.ipGlobal());
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

private:
    // material properties of the porous medium
    std::shared_ptr<SpatialParams> spatialParams_;
};

} // end namespace Dumux

#endif
