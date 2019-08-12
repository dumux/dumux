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
 * \ingroup FacetTests
 * \brief The problem for the lagrange-multiplier domain that is
 *        used to model contact mechanics of closing fractures.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_LAGRANGE_MP_PROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_LAGRANGE_MP_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dumux/common/feproblem.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/discretization/fem.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/multidomain/lagrangemultiplier/localresidual.hh>

namespace Dumux {

// The model traits corresponding to this problem
template<int dimWorld>
struct ContactLagrangeModelTraits
{
    static constexpr int numEq() { return dimWorld; }
};

// forward declarations
template<class TypeTag> class LagrangeProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct LagrangeFacet { using InheritsFrom = std::tuple<FiniteElementModel, ModelProperties>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LagrangeFacet> { using type = Dune::FoamGrid<1, 2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::LagrangeFacet> { using type = LagrangeProblem<TypeTag>; };

// We use a lagrange basis of zero-th order here
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::LagrangeFacet>
{
private:
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, 0>;
public:
    using type = FEGridGeometry<FEBasis>;
};

// We use a lagrange basis of zero-th order here
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::LagrangeFacet>
{ using type = GetPropType<TypeTag, Properties::GridGeometry>; };

// The model traits specifying the number of equations
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::LagrangeFacet>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
public:
    using type = ContactLagrangeModelTraits<GridView::dimensionworld>;
};

// The local residual type
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::LagrangeFacet>
{ using type = LagrangeMultiplierLocalResidual<TypeTag>; };

} // end namespace Properties

/*!
 * \ingroup FacetTests
 * \brief The problem for the (d-1)-dimensional facet domain in the elastic
 *        single-phase facet coupling test.
 */
template<class TypeTag>
class LagrangeProblem : public FEProblem<TypeTag>
{
    using ParentType = FEProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SecondaryVariables = typename GridVariables::SecondaryVariables;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int numEq = NumEqVector::size();
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    LagrangeProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManagerPtr,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName")
                         + "_"
                         + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        initialValues_ = getParamFromGroup<NumEqVector>(paramGroup, "InitialLagrangeMP");
        penaltyFactor_ = getParamFromGroup<Scalar>(paramGroup, "PenaltyFactor");
        frictionCoefficient_ = getParamFromGroup<Scalar>(paramGroup, "FrictionCoefficient");
        initialAperture_ = getParam<Scalar>("OnePFacet.SpatialParams.InitialAperture");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_; }

    //! Specifies the kind of boundary condition at a boundary position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // the tangential traction must be zero at fracture tips
        // values.setDirichlet(1);

        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     */
    template<class ElementSolution, class IpData>
    NumEqVector evalConstraint(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& elemSol,
                               const IpData& ipData,
                               const SecondaryVariables& secVars) const
    {
        using std::max;
        NumEqVector result;

        const auto& traction = secVars.priVars();
        const auto& contactSurfaceSegment = couplingManager().getContactSurfaceSegment(element);

        const auto a = couplingManager().computeAperture(element, ipData.ipGlobal(), initialAperture_);
        const auto deltaUT = couplingManager().computeTangentialDisplacementJump(element, ipData.ipGlobal());

        const auto tN = traction*contactSurfaceSegment.getBasisVector(dimWorld-1);
        auto tT = contactSurfaceSegment.getBasisVector(dimWorld-1);
        tT *= -1.0;
        tT *= tN;
        tT += traction;
        // if (feGeometry.gridGeometry().elementMapper().index(element) > -1){
        // std::cout << "ip: " << ipData.ipGlobal() << " -> a " << a << " -> delta " << deltaUT << " // "
        //           << couplingManager().computeNormalDisplacementJump(element, ipData.ipGlobal()) << std::endl;
        // std::cout << "TRACTION " << traction << std::endl;
        // std::cout << "BASIS Vector1:" << contactSurfaceSegment.getBasisVector(0) << std::endl;
        // std::cout << "BASIS Vector2:" << contactSurfaceSegment.getBasisVector(dimWorld-1) << std::endl;
        // }
        const auto normalMaxArg = -1.0*tN - penaltyFactor_*a;
        auto tangMaxArg = deltaUT;
        tangMaxArg *= penaltyFactor_;
        tangMaxArg -= tT;
        // if (feGeometry.gridGeometry().elementMapper().index(element) > -1){
        // std::cout << "NORMALMAXARG " << normalMaxArg << std::endl;
        // std::cout << "TANGMAXARG " << tangMaxArg << std::endl;
        // }

        const auto tangMaxArgNorm = tangMaxArg.two_norm();
        const auto& tangent = contactSurfaceSegment.getBasisVector(0);
        const auto fricCoeff = frictionCoefficientAtPos(ipData.ipGlobal());
        result[0] = -1.0*(tT*tangent)*max(fricCoeff*normalMaxArg, tangMaxArgNorm)
                    -fricCoeff*max(0.0, normalMaxArg)*(tangMaxArg*tangent);
        result[numEq-1] = -1.0*tN - max(0.0, normalMaxArg);
        // if (feGeometry.gridGeometry().elementMapper().index(element) > -1){
        // std::cout << "RESIDUAL = " << result << std::endl;
    // }
        return result;
    }

    //! Sets the aperture as extrusion factor.
    template<class IpData, class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const IpData& ipData,
                           const ElementSolution& elemSol) const
    { return 1.0; } //return couplingManager().computeAperture(element, scv, initialAperture_); }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return initialValues_; }

    //! Evaluates the Dirichlet boundary conditions.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables({0.0, 0.0}); }

    //! Returns const reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! Returns the friction coefficient at a given position
    const Scalar frictionCoefficientAtPos(const GlobalPosition& globalPos) const
    {
        auto f = frictionCoefficient_;
        double d;
        if (globalPos[1] > 0.5)
            d = (GlobalPosition({0.5, 0.75}) - globalPos).two_norm2();
        else
            d = (GlobalPosition({0.5, 0.25}) - globalPos).two_norm2();

        static const auto factor = getParam<Scalar>("Problem.FrictionCoefficientParam");
        return f*(1.0 + factor*std::exp(-800.0*d));
    }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::string problemName_;

    PrimaryVariables initialValues_;
    Scalar penaltyFactor_;
    Scalar frictionCoefficient_;
    Scalar initialAperture_;
};

} // end namespace Dumux

#endif
