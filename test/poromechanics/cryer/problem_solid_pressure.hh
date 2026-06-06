// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Solid pressure subdomain problem for the Cryer sphere benchmark.
 *
 * Enforces the algebraic volumetric constitutive law:
 *   p_s^eq(J) = K * (J^2 - 1) / (2J)
 *
 * where K is the drained bulk modulus.  For K = 1 MPa (nu=0 material), this gives
 * a non-zero volumetric solid stress that, together with the deviatoric Neo-Hookean
 * stress, reproduces the correct linearised Biot elastic response.
 */
#ifndef DUMUX_CRYER_PROBLEM_SOLID_PRESSURE_HH
#define DUMUX_CRYER_PROBLEM_SOLID_PRESSURE_HH

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

template<class TypeTag>
class CryerSolidPressureProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryFace = typename GridGeometry::BoundaryFace;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    CryerSolidPressureProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                               std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    /*!
     * \brief Volumetric pressure from the strain energy density.
     *
     * For ψ_vol = K/4*(J^2-1) - K/2*lnJ (Neo-Hookean volumetric-isochoric split):
     *   p_s^eq = dψ_vol/dJ = K*(J^2-1)/(2J)
     *
     * At J=1: p_s^eq = 0 (stress-free reference configuration). ✓
     *
     * Note: In the Neo-Hookean VI model, the parameter in W_vol is the drained bulk
     * modulus K directly (not the standard Lamé lambda = K - 2G/3).  The linearised
     * drained stress is σ = 2G·E_dev + K·ε_v·I (same as standard linear elasticity
     * with the same K and G), because P_iso → 2G·E_dev (zero volumetric part) and
     * P_vol → K·ε_v·I.
     */
    Scalar volumetricPressure(Scalar J) const
    {
        const Scalar K = this->spatialParams().bulkModulus();
        return K * (J*J - 1.0) / (2.0*J);
    }

    //! New boundary types interface — solid pressure has no Dirichlet BCs
    BoundaryTypes boundaryTypes(const auto& /*fvGeometry*/, const BoundaryFace&) const
    {
        BoundaryTypes bct;
        bct.setAllFluxBoundary();   // purely algebraic constraint, never Dirichlet
        return bct;
    }

    template<class ElementVariables>
    NumEqVector boundaryFluxIntegral(const auto&, const ElementVariables&, const auto&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux
#endif
