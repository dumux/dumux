// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Flow of a slightly compressible fluid in the pore network of a deforming sphere packing (requires UMFPack)
 */
#ifndef DUMUX_PNM_EXTRACTION_PORE_SCALE_FLOW_HH
#define DUMUX_PNM_EXTRACTION_PORE_SCALE_FLOW_HH

#include <array>
#include <memory>
#include <optional>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/porenetwork/extraction/spherepackingnetwork.hh>

namespace Dumux::PoreNetwork::SpherePacking {

/*!
 * \brief Pore pressures of a fluid in the network of a deforming sphere packing
 *
 * In every pore without imposed pressure the outflow balances the volume change of the pore, a volume
 * source and, for a fluid of finite bulk modulus K_f, the compression of the fluid in the pore:
 * (V_i/K_f) dp_i/dt + sum_j g_ij (p_i - p_j) = -dV_i/dt + q_i, with the void volume V_i, the conductances
 * g_ij = A R_h^2/(2 mu L) of the throats and the volume source q_i (pore-scale finite volume method,
 * Chareyre et al. 2012, Catalano et al. 2014). The storage term is discretised implicitly over the time
 * step; without bulk modulus the fluid is incompressible. Pores touching a wall with an imposed pressure
 * take that pressure. The matrix depends on the network and the time step only; it is factorised once,
 * and each solve is a back substitution.
 */
template<class Scalar>
class PoreScaleFlow
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

public:
    /*!
     * \param network the network of a packing with walls (pore labels are wall masks)
     * \param viscosity dynamic viscosity of the fluid
     * \param wallPressure the imposed pressure of each wall, none for an impermeable wall
     * \param bulkModulus bulk modulus of the fluid, none for an incompressible fluid
     * \param dt time step of the storage term
     */
    PoreScaleFlow(const Network<Scalar>& network, Scalar viscosity,
                  const std::array<std::optional<Scalar>, 6>& wallPressure,
                  std::optional<Scalar> bulkModulus = std::nullopt, Scalar dt = 1.0)
    : network_(network)
    , dirichlet_(network.pores.size())
    , storage_(network.pores.size(), 0.0)
    {
        if (bulkModulus)
            for (std::size_t i = 0; i < network.pores.size(); ++i)
                storage_[i] = network.pores[i].volume/(*bulkModulus*dt);

        const std::size_t n = network.pores.size();
        for (std::size_t i = 0; i < n; ++i)
            for (int k = 0; k < 6; ++k)
                if (wallPressure[k] && network.pores[i].label > 0 && (network.pores[i].label & (1 << k)))
                    dirichlet_[i] = *wallPressure[k];

        Dune::MatrixIndexSet pattern(n, n);
        for (std::size_t i = 0; i < n; ++i)
            pattern.add(i, i);
        for (const auto& throat : network.throats)
        {
            pattern.add(throat.pores[0], throat.pores[1]);
            pattern.add(throat.pores[1], throat.pores[0]);
        }
        matrix_ = std::make_unique<Matrix>();
        pattern.exportIdx(*matrix_);
        *matrix_ = 0.0;

        for (const auto& throat : network.throats)
        {
            const Scalar g = TransmissibilityChareyre<Scalar>::singlePhaseTransmissibility(
                throat.fluidArea, throat.hydraulicRadius, throat.length)/viscosity;
            const auto i = throat.pores[0], j = throat.pores[1];
            (*matrix_)[i][i] += g;
            (*matrix_)[j][j] += g;
            (*matrix_)[i][j] -= g;
            (*matrix_)[j][i] -= g;
        }
        for (std::size_t i = 0; i < n; ++i)
            (*matrix_)[i][i] += storage_[i];
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!dirichlet_[i])
                continue;
            for (auto entry = (*matrix_)[i].begin(); entry != (*matrix_)[i].end(); ++entry)
                *entry = entry.index() == i ? 1.0 : 0.0;
        }

        solver_ = std::make_unique<Dune::UMFPack<Matrix>>(*matrix_, 0);
        pressure_.resize(n, 0.0);
    }

    /*!
     * \brief Pore pressures at the end of a time step for the given rates of change of the pore volumes and
     *        volume sources, starting from the current pressures
     */
    const std::vector<Scalar>& solve(const std::vector<Scalar>& volumeRates, const std::vector<Scalar>& sources = {})
    {
        const std::size_t n = network_.pores.size();
        Vector rhs(n), x(n);
        for (std::size_t i = 0; i < n; ++i)
            rhs[i] = dirichlet_[i] ? *dirichlet_[i]
                                   : -volumeRates[i] + (sources.empty() ? 0.0 : sources[i]) + storage_[i]*pressure_[i];
        Dune::InverseOperatorResult result;
        solver_->apply(x, rhs, result);
        for (std::size_t i = 0; i < n; ++i)
            pressure_[i] = x[i];
        return pressure_;
    }

    const std::vector<Scalar>& pressure() const
    { return pressure_; }

    //! pressures to start the next time step from, e.g. transferred from a previous network
    void setPressure(const std::vector<Scalar>& pressure)
    { pressure_ = pressure; }

    //! whether the pressure of a pore is imposed
    bool isDirichlet(std::size_t pore) const
    { return dirichlet_[pore].has_value(); }

private:
    const Network<Scalar>& network_;
    std::vector<std::optional<Scalar>> dirichlet_;
    std::vector<Scalar> storage_;
    std::unique_ptr<Matrix> matrix_;
    std::unique_ptr<Dune::UMFPack<Matrix>> solver_;
    std::vector<Scalar> pressure_;
};

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
