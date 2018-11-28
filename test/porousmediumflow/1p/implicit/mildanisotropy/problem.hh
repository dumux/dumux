// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup OnePTests
 * \brief The problem for the incompressible single-phase
 *        flow test based on "test 1: mild anisotropy" from:
 *
 *        Benchmark on Discretization Schemes for
 *        Anisotropic Diffusion Problems on General Grids
 *        (https://hal.archives-ouvertes.fr/hal-00429843/)
 */
#ifndef DUMUX_ONEP_MILDANISOTROPY_TEST_PROBLEM_HH
#define DUMUX_ONEP_MILDANISOTROPY_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class MildAnisotropyProblem;

namespace Properties {
// create the type tag nodes
// Create new type tags
namespace TTag {
struct MildAnisotropy { using InheritsFrom = std::tuple<OneP>; };
struct MildAnisotropyTpfa { using InheritsFrom = std::tuple<MildAnisotropy, CCTpfaModel>; };
struct MildAnisotropyMpfa { using InheritsFrom = std::tuple<MildAnisotropy, CCMpfaModel>; };
struct MildAnisotropyBox { using InheritsFrom = std::tuple<MildAnisotropy, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MildAnisotropy> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::MildAnisotropy> { using type = MildAnisotropyProblem<TypeTag>; };

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::MildAnisotropy> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MildAnisotropy>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = MildAnisotropySpatialParams<FVGridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MildAnisotropy>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief The problem for the incompressible single-phase
 *        flow test based on "test 1: mild anisotropy" from:
 *
 *        Benchmark on Discretization Schemes for
 *        Anisotropic Diffusion Problems on General Grids
 *        (https://hal.archives-ouvertes.fr/hal-00429843/)
 */
template<class TypeTag>
class MildAnisotropyProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

public:
    MildAnisotropyProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    , considerNeumannBC_(getParam<bool>("Problem.ConsiderNeumannBC"))
    {}

    //! Returns the exact solution at a given position
    Scalar exactSolution(const GlobalPosition& globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return 16.0*x*(1.0-x)*y*(1.0-y);
    }

    //! Returns the exact solution gradient at a given position
    GlobalPosition exactGradient(const GlobalPosition& globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        GlobalPosition grad;
        grad[0] = y*(1.0-y)*(16.0-32.0*x);
        grad[1] = 16.0*x*(1.0-x)*(1.0-2.0*y);
        return grad;
    }

    //! Specify the type of boundary conditions for a position on the boundary
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        if (considerNeumannBC_ && (globalPos[0] < 1e-6 ||
                                   globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - 1e-6))
            values.setAllNeumann();
        return values;
    }

    //! evaluate Dirichlet boundary conditions at a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables( exactSolution(globalPos) ); }

    //! evaluate the source term in a control volume
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        const auto& ip = scv.center();
        const auto x = ip[0];
        const auto y = ip[1];
        const Scalar dgradUx_dx = -32*y*(1.0-y);
        const Scalar dgradUy_dx = (16.0-32.0*x)*(1.0-2*y);
        const Scalar dgradUx_dy = dgradUy_dx;
        const Scalar dgradUy_dy = -32*x*(1.0-x);

        Scalar source = 0.0;
        const auto K = elemVolVars[scv].permeability();
        source -= K[0][0]*dgradUx_dx + K[0][1]*dgradUy_dx;
        source -= K[1][0]*dgradUx_dy + K[1][1]*dgradUy_dy;

        return NumEqVector(source);
    }

    //! evaluate Neumann boundary conditions on a given boundary segment
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto K = elemVolVars[insideScv].permeability();

        GlobalPosition KgradU;
        K.mv(exactGradient(scvf.ipGlobal()), KgradU);

        auto flux = KgradU*scvf.unitOuterNormal();
        flux *= -1.0;

        return NumEqVector(flux);
    }

    //! returns the temperature in the domain
    Scalar temperature() const
    { return 283.15; /* 10°C*/ }

private:
    bool considerNeumannBC_;
};

} // end namespace Dumux

#endif
