// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The problem for the bulk domain in the single-phase facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULKPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULKPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief Test problem for the incompressible one-phase model
 *        with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    OnePBulkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManagerPtr,
                    const std::string& paramGroup = "Bulk")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    , lowDimPermeability_(getParam<Scalar>("LowDim.SpatialParams.Permeability"))
    , aperture_(getParam<Scalar>("LowDim.SpatialParams.Aperture"))
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::cosh;
        Scalar u = (1.0 - lowDimPermeability_)*cos(globalPos[0])*cosh(aperture_/2);
        return NumEqVector(u);
    }

    //! Evaluates the exact solution at a given position.
    Scalar exact(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::cosh;
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return lowDimPermeability_*cos(x)*cosh(y) + (1.0 - lowDimPermeability_)*cos(x)*cosh(aperture_/2);
    }

    //! Evaluates the exact gradient at a given position.
    GlobalPosition exactGradient(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::sin;
        using std::cosh;
        using std::sinh;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        GlobalPosition gradU;
        gradU[0] = -lowDimPermeability_*sin(x)*cosh(y) + (lowDimPermeability_ - 1.0)*sin(x)*cosh(aperture_/2);
        gradU[1] = lowDimPermeability_*cos(x)*sinh(y);

        return gradU;
    }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(exact(globalPos)); }

    //! Evaluates the Neumann boundary condition for a boundary segment.
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        auto pos = scvf.ipGlobal();
        const Scalar k = this->spatialParams().permeabilityAtPos(pos);
        const auto gradU = exactGradient(pos);
        return NumEqVector( -1.0*k*(gradU*scvf.unitOuterNormal()) );
    }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar lowDimPermeability_;
    Scalar aperture_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
