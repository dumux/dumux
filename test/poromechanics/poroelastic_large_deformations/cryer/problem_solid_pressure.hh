// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Solid pressure subdomain problem for the Cryer sphere benchmark.
 */
#ifndef DUMUX_TEST_POROELASTIC_LARGE_DEF_CRYER_PROBLEM_SOLID_PRESSURE_HH
#define DUMUX_TEST_POROELASTIC_LARGE_DEF_CRYER_PROBLEM_SOLID_PRESSURE_HH

#include <memory>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

template<class TypeTag>
class CryerSolidPressureProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryFace = typename GridDiscretization::BoundaryFace;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    CryerSolidPressureProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                               std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    /*!
     * \brief Volumetric pressure \f$ p_s = \partial\psi_\text{vol}/\partial J = K (J^2 - 1)/(2J) \f$
     *        of \f$ \psi_\text{vol} = \frac{K}{4}(J^2 - 1) - \frac{K}{2}\ln J \f$
     *
     * The linearization \f$ p_s \approx K (J - 1) \f$ makes \f$ K \f$ the drained bulk modulus
     * of the linear theory of Cryer's solution.
     */
    Scalar volumetricPressure(Scalar J) const
    {
        const Scalar K = this->spatialParams().bulkModulus();
        return K * (J*J - 1.0) / (2.0*J);
    }

    //! The solid pressure equation is algebraic and has no boundary conditions
    BoundaryTypes boundaryTypes(const auto& /*elemDisc*/, const BoundaryFace&) const
    {
        BoundaryTypes bct;
        bct.setAllFluxBoundary();
        return bct;
    }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux
#endif
