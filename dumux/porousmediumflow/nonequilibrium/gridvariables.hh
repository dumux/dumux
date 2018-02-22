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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief Class storing scv and scvf variables
 */
#ifndef DUMUX_NONEQUILIBRIUM_GRID_VARIABLES_HH
#define DUMUX_NONEQUILIBRIUM_GRID_VARIABLES_HH

#include <dune/common/fvector.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvgridvariables.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief This class stores the velocities which are used to compute reynoldsnumbers for the source terms of nonequilibrium models
 */
template<class TypeTag>
class NonEquilibriumGridVariables: public FVGridVariables<TypeTag>
{
    using ParentType = FVGridVariables<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VelocityOutput = typename GET_PROP_TYPE(TypeTag, VelocityOutput);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum { dim = GridView::dimension }; // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };

    using GlobalPosition = Dune::FieldVector<typename GridView::Grid::ctype, dimWorld>;
    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethod::box;

public:
    //! Constructor
    NonEquilibriumGridVariables(std::shared_ptr<const Problem> problem,
                                std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(problem, fvGridGeometry)
    , problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    {
        for (int phaseIdx =0; phaseIdx<numPhases; ++phaseIdx)
        {
            velocityNorm_[phaseIdx].assign(fvGridGeometry->numDofs(), 0.0);
        }
    }

    void calcVelocityAverage(const SolutionVector& curSol)
    {
        // instatiate the velocity output
        VelocityOutput velocityOutput(*problem_, *fvGridGeometry_, *this, curSol);
        std::array<std::vector<GlobalPosition>, numPhases> velocity;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            if(isBox && dim == 1)
                velocity[phaseIdx].resize(fvGridGeometry_->gridView().size(0));
            else
                velocity[phaseIdx].resize(fvGridGeometry_->numDofs());
        }
        for (const auto& element : elements(fvGridGeometry_->gridView(), Dune::Partitions::interior))
        {
            const auto eIdxGlobal = fvGridGeometry_->elementMapper().index(element);

            auto fvGeometry = localView(*fvGridGeometry_);
            auto elemVolVars = localView(this->curGridVolVars());

            fvGeometry.bind(element);

            elemVolVars.bind(element, fvGeometry, curSol);

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                velocityOutput.calculateVelocity(velocity[phaseIdx], elemVolVars, fvGeometry, element, phaseIdx);

                for (auto&& scv : scvs(fvGeometry))
                {
                    const auto dofIdxGlobal = scv.dofIndex();
                    if (isBox && dim == 1)
                        velocityNorm_[phaseIdx][dofIdxGlobal] = velocity[phaseIdx][eIdxGlobal].two_norm();
                    else
                       velocityNorm_[phaseIdx][dofIdxGlobal] = velocity[phaseIdx][dofIdxGlobal].two_norm();
                }
            }//end phases
        }//end elements
  }// end calcVelocity

    /*!
     * \brief Access to the averaged (magnitude of) velocity for each vertex.
     *
     * \param phaseIdx The index of the fluid phase
     * \param dofIdxGlobal The global index of the degree of freedom
     *
     */
    const Scalar volumeDarcyMagVelocity(const unsigned int phaseIdx,
                                        const unsigned int dofIdxGlobal) const
    { return velocityNorm_[phaseIdx][dofIdxGlobal]; }

private:
    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;
    std::array<std::vector<Dune::FieldVector<Scalar, 1> > , numPhases> velocityNorm_;


};
} // end namespace

#endif
