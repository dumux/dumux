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
 * \ingroup TracerTests
 * \brief The properties for the incompressible test
 */
#ifndef DUMUX_ONEP_TRACER_TEST_PROBLEM_HH
#define DUMUX_ONEP_TRACER_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "spatialparams_1p.hh"

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief The properties for the incompressible test
 */
// forward declarations
template<class TypeTag>
class OnePTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct IncompressibleTest { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::IncompressibleTest> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::IncompressibleTest> { using type = OnePTestProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::IncompressibleTest>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::IncompressibleTest> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IncompressibleTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
} // end namespace Properties

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        const auto globalPos = scvf.ipGlobal();

        Scalar eps = 1.0e-6;
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        const auto& pos = scvf.ipGlobal();
        PrimaryVariables values(0);
        values[0] = 1.0e+5*(1.1 - pos[dimWorld-1]*0.1);
        return values;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

};

} // end namespace Dumux

#endif
