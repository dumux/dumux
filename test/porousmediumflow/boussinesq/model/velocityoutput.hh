// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Velocity output for the vector-potential Boussinesq model (2D/3D).
 *
 * Recovers the Darcy velocity u = ∇ × A by averaging the curl of the
 * vector-potential gradient over all interior SCVFs of each element,
 * then distributes the element velocity to its vertices (inverse-count averaging).
 *
 * In 2D (single potential ψ):  u = (∂ψ/∂y, −∂ψ/∂x)
 * In 3D (three potentials A):  u = ∇ × A  (full curl)
 */
#ifndef DUMUX_BOUSSINESQ_VORTICITY_VELOCITY_OUTPUT_HH
#define DUMUX_BOUSSINESQ_VORTICITY_VELOCITY_OUTPUT_HH

#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/velocityoutput.hh>

namespace Dumux {

template<class GridVariables>
class BoussinesqVorticityVelocityOutput
    : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using GridGeometry          = typename GridVariables::GridGeometry;
    using GridView              = typename GridGeometry::GridView;
    using GridVolumeVariables   = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVarsCache  = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry     = typename GridGeometry::LocalView;
    using Element               = typename GridView::template Codim<0>::Entity;
    using Scalar                = typename GridVariables::Scalar;

    static constexpr int dim      = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Velocity = Dune::FieldVector<Scalar, dimWorld>;

    // Number of vector-potential components resolved in this dimension
    static constexpr int nPot = (dimWorld == 2) ? 1 : dimWorld;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    explicit BoussinesqVorticityVelocityOutput(const GridVariables& gridVariables)
    {
        enableOutput_ = getParamFromGroup<bool>(
            gridVariables.curGridVolVars().problem().paramGroup(),
            "Vtk.AddVelocity", false);

        if (enableOutput_)
        {
            const auto& gg = gridVariables.gridGeometry();
            cellNum_.assign(gg.gridView().size(dim), 0);
            for (const auto& element : elements(gg.gridView()))
                for (unsigned int vIdx = 0; vIdx < element.subEntities(dim); ++vIdx)
                    ++cellNum_[gg.vertexMapper().subIndex(element, vIdx, dim)];
        }
    }

    bool        enableOutput()   const override { return enableOutput_; }
    std::string phaseName(int)   const override { return "fluid"; }
    int         numFluidPhases() const override { return 1; }

    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int /*phaseIdx*/) const override
    {
        if (!enableOutput_) return;

        // Average the vector-potential gradients over all interior SCVFs
        std::array<Velocity, nPot> gradA;
        for (auto& g : gradA) g = Velocity(0.0);
        int nInteriorFaces = 0;

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary()) continue;
            ++nInteriorFaces;

            const auto& fc = elemFluxVarsCache[scvf];
            std::array<Velocity, nPot> localGrad;
            for (auto& g : localGrad) g = Velocity(0.0);

            for (const auto& scv : scvs(fvGeometry))
                for (int k = 0; k < nPot; ++k)
                    localGrad[k].axpy(elemVolVars[scv].vectorPotential(k),
                                      fc.gradN(scv.indexInElement()));

            for (int k = 0; k < nPot; ++k)
                gradA[k] += localGrad[k];
        }

        if (nInteriorFaces == 0) return;
        for (auto& g : gradA) g /= Scalar(nInteriorFaces);

        // Compute u = ∇ × A
        Velocity u(0.0);
        if constexpr (dimWorld == 2)
        {
            // ψ ≡ A[0]: u_x = ∂ψ/∂y,  u_y = -∂ψ/∂x
            u[0] =  gradA[0][1];
            u[1] = -gradA[0][0];
        }
        else
        {
            // A[0]=A_x, A[1]=A_y, A[2]=A_z
            u[0] = gradA[2][1] - gradA[1][2];
            u[1] = gradA[0][2] - gradA[2][0];
            u[2] = gradA[1][0] - gradA[0][1];
        }

        // Distribute to vertices weighted by inverse element-per-vertex count
        for (const auto& scv : scvs(fvGeometry))
        {
            const int vIdx = scv.dofIndex();
            velocity[vIdx] += u / Scalar(cellNum_[vIdx]);
        }
    }

private:
    bool             enableOutput_ = false;
    std::vector<int> cellNum_;
};

} // namespace Dumux

#endif
