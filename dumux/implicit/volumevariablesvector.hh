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
 * \brief Base class for the volume variables vector
 */
#ifndef DUMUX_IMPLICIT_VOLVARSVECTOR_HH
#define DUMUX_IMPLICIT_VOLVARSVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag>
class VolumeVariablesVector : public std::vector<typename GET_PROP_TYPE(TypeTag, VolumeVariables)>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using IndexType = typename GridView::IndexSet::IndexType;

    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    void update(const Problem& problem, const SolutionVector& sol)
    {
        numScvs_ = problem.model().fvGeometries().numScv();
        if (!isBox)
            numBoundaryScvs_ = problem.model().fvGeometries().numBoundaryScvf();
        else
            numBoundaryScvs_ = 0;

        volumeVariables_.resize(numScvs_);
        boundaryVolumeVariables_.resize(numBoundaryScvs_);
        for (const auto& element : elements(problem.gridView()))
        {
            for (auto&& scv : problem.model().fvGeometries(element).scvs())
            {
                (*this)[scv.index()].update(sol[scv.dofIndex()],
                                            problem,
                                            element,
                                            scv);
            }

            if (!isBox)
            {
                for (auto&& scvFace : problem.model().fvGeometries(element).scvfs())
                {
                    // if we are not on a boundary, skip the rest
                    if (!scvFace.boundary())
                        continue;

                    // When complex boundary handling is inactive, we only use BC vol vars on pure Dirichlet boundaries
                    auto bcTypes = problem.boundaryTypes(element, scvFace);
                    if (/*TODO !GET_PROP_VALUE(TypeTag, BoundaryReconstruction) && */!(bcTypes.hasDirichlet() && !bcTypes.hasNeumann()))
                        continue;

                    const auto insideScvIdx = scvFace.insideScvIdx();
                    const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
                    const auto dirichletPriVars = problem.dirichlet(element, scvFace);

                    (*this)[scvFace.outsideScvIdx()].update(dirichletPriVars, problem, element, insideScv);
                }
            }
        }
    }

    const VolumeVariables& operator [](IndexType scvIdx) const
    {
        assert(scvIdx < numScvs_ + numBoundaryScvs_);

        if (scvIdx < numScvs_)
            return volumeVariables_[scvIdx];
        else
            return boundaryVolumeVariables_[scvIdx-numScvs_];
    }

    VolumeVariables& operator [](IndexType scvIdx)
    {
        assert(scvIdx < numScvs_ + numBoundaryScvs_);

        if (scvIdx < numScvs_)
            return volumeVariables_[scvIdx];
        else
            return boundaryVolumeVariables_[scvIdx-numScvs_];
    }

private:
    IndexType numScvs_;
    IndexType numBoundaryScvs_;
    std::vector<VolumeVariables> volumeVariables_;
    std::vector<VolumeVariables> boundaryVolumeVariables_;
};

} // end namespace

#endif
