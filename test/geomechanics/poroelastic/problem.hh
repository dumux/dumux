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
/**
 * \file
 * \brief Definition of a test problem for the poro-elastic model
 */
#ifndef DUMUX_POROELASTIC_PROBLEM_HH
#define DUMUX_POROELASTIC_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class PoroElasticProblem;

namespace Properties {
NEW_TYPE_TAG(TestPoroElastic, INHERITS_FROM(BoxModel, PoroElastic));
// Set the grid type
SET_TYPE_PROP(TestPoroElastic, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(TestPoroElastic, Problem, Dumux::PoroElasticProblem<TypeTag>);
// The fluid phase consists of one constant component
SET_TYPE_PROP(TestPoroElastic,
              FluidSystem,
              Dumux::FluidSystems::OnePLiquid< typename GET_PROP_TYPE(TypeTag, Scalar),
                                               Dumux::Components::Constant<0, typename GET_PROP_TYPE(TypeTag, Scalar)> >);
// The spatial parameters property
SET_TYPE_PROP(TestPoroElastic, SpatialParams, PoroElasticSpatialParams< typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                           typename GET_PROP_TYPE(TypeTag, FVGridGeometry) >);
}

/*!
 * \ingroup Geomechanics
 * \ingroup PoroElastic
 *
 * \brief Problem definition for the deformation of a poro-elastic body
 */
template<class TypeTag>
class PoroElasticProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr Scalar pi = M_PI;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GradU = Dune::FieldMatrix<Scalar, dim, dimWorld>;

public:
    //! The constructor
    PoroElasticProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {}

    //! The temperature in the domain
    static constexpr Scalar temperature()
    { return 273.15; }

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluate the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the effective fluid density
     * \param globalPos The global position
     */
    Scalar effectiveFluidDensityAtPos(const GlobalPosition& globalPos) const
    {
        // This test uses the constant component, obtain density only once
        using FS = typename GET_PROP_TYPE(TypeTag, FluidSystem);
        static const Scalar rho = FS::density( effectivePorePressureAtPos(globalPos), temperature() );
        return rho;
    }

    /*!
     * \brief Returns the effective pore pressure
     * \note We use the x-displacement as pressure solution. The shift to
     *       higher values is done to see a mor pronounced effect in stresses.
     *
     * \param globalPos The global position
     */
    Scalar effectivePorePressureAtPos(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos)[0] + 10; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     */
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        using std::sin;
        using std::cos;

        const auto ipGlobal = scv.center();
        const auto x = ipGlobal[0];
        const auto y = ipGlobal[1];

        // the lame parameters (we know they only depend on position here)
        const auto& lameParams = this->spatialParams().lameParamsAtPos(scv.center());
        const auto lambda = lameParams.lambda();
        const auto mu = lameParams.mu();

        // precalculated products
        const Scalar pi_2 = 2.0*pi;
        const Scalar pi_2_square = pi_2*pi_2;
        const Scalar cos_2pix = cos(pi_2*x);
        const Scalar sin_2pix = sin(pi_2*x);
        const Scalar cos_2piy = cos(pi_2*y);
        const Scalar sin_2piy = sin(pi_2*y);

        const Scalar dE11_dx = -2.0*sin_2piy;
        const Scalar dE22_dx = pi_2_square*cos_2pix*cos_2piy;
        const Scalar dE11_dy = pi_2*(1.0-2.0*x)*cos_2piy;
        const Scalar dE22_dy = -1.0*pi_2_square*sin_2pix*sin_2piy;
        const Scalar dE12_dy = 0.5*pi_2_square*(cos_2pix*cos_2piy - (x-x*x)*sin_2piy);
        const Scalar dE21_dx = 0.5*((1.0-2*x)*pi_2*cos_2piy - pi_2_square*sin_2pix*sin_2piy);

        // compute exact divergence of sigma
        PrimaryVariables divSigma(0.0);
        divSigma[Indices::momentum(/*x-dir*/0)] = lambda*(dE11_dx + dE22_dx) + 2*mu*(dE11_dx + dE12_dy);
        divSigma[Indices::momentum(/*y-dir*/1)] = lambda*(dE11_dy + dE22_dy) + 2*mu*(dE21_dx + dE22_dy);
        return divSigma;
    }

    /*!
     * \brief Evaluate the exact displacement to this problem at a given position.
     */
    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    {
        using std::sin;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        PrimaryVariables exact(0.0);
        exact[Indices::momentum(/*x-dir*/0)] = (x-x*x)*sin(2*pi*y);
        exact[Indices::momentum(/*y-dir*/1)] = sin(2*pi*x)*sin(2*pi*y);
        return exact;
    }

    /*!
     * \brief Evaluate the exact displacement gradient to this problem at a given position.
     */
    GradU exactGradient(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        static constexpr int xIdx = Indices::momentum(/*x-dir*/0);
        static constexpr int yIdx = Indices::momentum(/*y-dir*/1);

        GradU exactGrad(0.0);
        exactGrad[xIdx][xIdx] = (1-2*x)*sin(2*pi*y);
        exactGrad[xIdx][yIdx] = (x - x*x)*2*pi*cos(2*pi*y);
        exactGrad[yIdx][xIdx] = 2*pi*cos(2*pi*x)*sin(2*pi*y);
        exactGrad[yIdx][yIdx] = 2*pi*sin(2*pi*x)*cos(2*pi*y);
        return exactGrad;
    }

private:
    static constexpr Scalar eps_ = 3e-6;
};

} //end namespace

#endif
