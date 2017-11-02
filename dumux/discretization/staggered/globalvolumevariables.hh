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
 * \brief The global volume variables class for cell centered models
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GLOBAL_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GLOBAL_VOLUMEVARIABLES_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/discretization/staggered/elementvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class StaggeredGlobalVolumeVariables
{};

//! specialization in case of storing the volume variables
template<class TypeTag>
class StaggeredGlobalVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/true>
{
    // The local class needs to access and change volVars
    friend StaggeredElementVolumeVariables<TypeTag, true>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    enum { numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter) };

public:
    StaggeredGlobalVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        auto numScv = fvGridGeometry.numScv();
        auto numBoundaryScvf = fvGridGeometry.numBoundaryScvf();

        volumeVariables_.resize(numScv + numBoundaryScvf);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                PrimaryVariables priVars(0.0);
                priVars[cellCenterIdx] = sol[cellCenterIdx][scv.dofIndex()];
                ElementSolutionVector elemSol{std::move(priVars)};
                volumeVariables_[scv.dofIndex()].update(elemSol, problem(), element, scv);
            }

            // handle the boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip the rest
                if (!scvf.boundary())
                    continue;

                const auto bcTypes = problem().boundaryTypes(element, scvf);
                const auto insideScvIdx = scvf.insideScvIdx();
                const auto& insideScv = fvGeometry.scv(insideScvIdx);

                PrimaryVariables boundaryPriVars(0.0);

                for(int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
                {
                    if(bcTypes.isDirichlet(eqIdx) || bcTypes.isDirichletCell(eqIdx))
                        boundaryPriVars[cellCenterIdx][eqIdx] = problem().dirichlet(element, scvf)[cellCenterIdx][eqIdx];
                    else if(bcTypes.isNeumann(eqIdx) || bcTypes.isOutflow(eqIdx) || bcTypes.isSymmetry())
                        boundaryPriVars[cellCenterIdx][eqIdx] = sol[cellCenterIdx][scvf.insideScvIdx()][eqIdx];
                    //TODO: this assumes a zero-gradient for e.g. the pressure on the boundary
                    // could be made more general by allowing a non-zero-gradient, provided in problem file
                    else
                        if(eqIdx == Indices::pressureIdx)
                            DUNE_THROW(Dune::InvalidStateException, "Face at: " << scvf.center() << " has neither Dirichlet nor Neumann BC.");
                }
                ElementSolutionVector elemSol{std::move(boundaryPriVars)};
                volumeVariables_[scvf.outsideScvIdx()].update(elemSol, problem(), element, insideScv);
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementVolumeVariables localView(const StaggeredGlobalVolumeVariables& global)
    { return ElementVolumeVariables(global); }

    const VolumeVariables& volVars(const IndexType scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const IndexType scvIdx)
    { return volumeVariables_[scvIdx]; }

    const VolumeVariables& volVars(const SubControlVolume scv) const
    { return volumeVariables_[scv.dofIndex()]; }

    VolumeVariables& volVars(const SubControlVolume scv)
    { return volumeVariables_[scv.dofIndex()]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:

    const Problem* problemPtr_;

    std::vector<VolumeVariables> volumeVariables_;
};


//! Specialization when the current volume variables are not stored globally
template<class TypeTag>
class StaggeredGlobalVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/false>
{
    // local class needs access to the problem
    friend StaggeredElementVolumeVariables<TypeTag, false>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    void update(Problem& problem, const SolutionVector& sol)
    { problemPtr_ = &problem; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementVolumeVariables localView(const StaggeredGlobalVolumeVariables& global)
    { return ElementVolumeVariables(global); }

private:
    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;
};

} // end namespace

#endif
