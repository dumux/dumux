// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NonEquilibriumModel
 * \brief Class storing scv and scvf variables.
 */

#ifndef DUMUX_NONEQUILIBRIUM_GRID_VARIABLES_HH
#define DUMUX_NONEQUILIBRIUM_GRID_VARIABLES_HH

#include <memory>
#include <dune/common/fvector.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvgridvariables.hh>
#include <dumux/porousmediumflow/velocity.hh>
#include <dumux/porousmediumflow/fluxvariables.hh>

namespace Dumux {

/*!
 * \ingroup NonEquilibriumModel
 * \brief This class stores the velocities which are used to compute Reynolds
 *        numbers for the source terms of non-equilibrium models.
 */
template<class TypeTag>
class NonEquilibriumGridVariables
: public FVGridVariables<GetPropType<TypeTag, Properties::GridGeometry>,
                         GetPropType<TypeTag, Properties::GridVolumeVariables>,
                         GetPropType<TypeTag, Properties::GridFluxVariablesCache>>
{
    using ThisType = NonEquilibriumGridVariables<TypeTag>;
    using ParentType = FVGridVariables<GetPropType<TypeTag, Properties::GridGeometry>,
                                       GetPropType<TypeTag, Properties::GridVolumeVariables>,
                                       GetPropType<TypeTag, Properties::GridFluxVariablesCache>>;

    using VelocityBackend = PorousMediumFlowVelocity<ThisType, PorousMediumFluxVariables<TypeTag>>;

    static constexpr auto dim = ParentType::GridGeometry::GridView::dimension; // Grid and world dimension
    static constexpr auto dimWorld = ParentType::GridGeometry::GridView::dimensionworld;
    static constexpr int numPhases = ParentType::VolumeVariables::numFluidPhases();
    static constexpr bool isBox = ParentType::GridGeometry::discMethod == DiscretizationMethods::box;

public:
    //! Export the type used for scalar values
    using typename ParentType::Scalar;
    using typename ParentType::GridGeometry;

    //! Constructor
    template<class Problem>
    NonEquilibriumGridVariables(std::shared_ptr<Problem> problem,
                                std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(problem, gridGeometry)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            velocityNorm_[phaseIdx].assign(gridGeometry->numDofs(), 0.0);

        velocityBackend_ = std::make_unique<VelocityBackend>(*this);
    }

    template<class SolutionVector>
    void calcVelocityAverage(const SolutionVector& curSol)
    {
        using Scalar = typename SolutionVector::field_type;
        using VelocityVector = typename Dune::FieldVector<Scalar, dimWorld>;

        std::array<std::vector<VelocityVector>, numPhases> velocity;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            if(isBox && dim == 1)
                velocity[phaseIdx].resize(this->gridGeometry_->gridView().size(0));
            else
                velocity[phaseIdx].resize(this->gridGeometry_->numDofs());
        }

        auto fvGeometry = localView(*this->gridGeometry_);
        auto elemVolVars = localView(this->curGridVolVars());
        auto elemFluxVarsCache = localView(this->gridFluxVarsCache());

        for (const auto& element : elements(this->gridGeometry_->gridView(), Dune::Partitions::interior))
        {
            const auto eIdxGlobal = this->gridGeometry_->elementMapper().index(element);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                velocityBackend_->calculateVelocity(velocity[phaseIdx], element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);

                for (auto&& scv : scvs(fvGeometry))
                {
                    const auto dofIdxGlobal = scv.dofIndex();
                    if (isBox && dim == 1)
                        velocityNorm_[phaseIdx][dofIdxGlobal] = velocity[phaseIdx][eIdxGlobal].two_norm();
                    else
                       velocityNorm_[phaseIdx][dofIdxGlobal] = velocity[phaseIdx][dofIdxGlobal].two_norm();
                }
            } //end phases
        } //end elements
    } // end calcVelocity

    /*!
     * \brief Access to the averaged (magnitude of) velocity for each vertex.
     *
     * \param phaseIdx The index of the fluid phase
     * \param dofIdxGlobal The global index of the degree of freedom
     */
    const Scalar volumeDarcyMagVelocity(const unsigned int phaseIdx,
                                        const unsigned int dofIdxGlobal) const
    { return velocityNorm_[phaseIdx][dofIdxGlobal]; }

private:
    std::array<std::vector<Dune::FieldVector<Scalar, 1> > , numPhases> velocityNorm_;
    std::unique_ptr<VelocityBackend> velocityBackend_;
};

} // end namespace Dumux

#endif
