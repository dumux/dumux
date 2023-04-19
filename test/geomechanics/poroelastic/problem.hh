// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Definition of a test problem for the poro-elastic model.
 */
#ifndef DUMUX_POROELASTIC_PROBLEM_HH
#define DUMUX_POROELASTIC_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/geomechanics/fvproblem.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Problem definition for the deformation of a poro-elastic body.
 */
template<class TypeTag>
class PoroElasticProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr Scalar pi = M_PI;

public:
    PoroElasticProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
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

        // The source term is composed of the divergence of the stress tensor
        // resulting from the exact solution minus the pressure gradient
        PrimaryVariables divSigma(0.0);
        divSigma[Indices::momentum(/*x-dir*/0)] = lambda*(dE11_dx + dE22_dx) + 2*mu*(dE11_dx + dE12_dy);
        divSigma[Indices::momentum(/*y-dir*/1)] = lambda*(dE11_dy + dE22_dy) + 2*mu*(dE21_dx + dE22_dy);
        divSigma -= this->spatialParams().effectivePorePressureGradient(ipGlobal);
        return divSigma;
    }

    /*!
     * \brief Evaluates the exact displacement to this problem at a given position.
     */
    PrimaryVariables exactDisplacement(const GlobalPosition& globalPos) const
    {
        using std::sin;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        PrimaryVariables exact(0.0);
        exact[Indices::momentum(/*x-dir*/0)] = (x-x*x)*sin(2*pi*y);
        exact[Indices::momentum(/*y-dir*/1)] = sin(2*pi*x)*sin(2*pi*y);
        return exact;
    }

private:
    static constexpr Scalar eps_ = 3e-6;
};

} // end namespace Dumux

#endif
