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
 * \brief The global volume variables class
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GLOBAL_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_BOX_GLOBAL_VOLUMEVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/box/elementvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class BoxGlobalVolumeVariables
{};

// specialization in case of storing the volume variables
template<class TypeTag>
class BoxGlobalVolumeVariables<TypeTag,/*enableGlobalVolVarCache*/true>
{
    // The local class needs to access and change volVars
    friend BoxElementVolumeVariables<TypeTag, true>;
    // The local jacobian needs to access and change volVars for derivative calculation
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;

        volumeVariables_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem_().elementMapper().index(element);

            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bindElement(element);

            // get the element solution
            auto elemSol = problem.model().elementSolution(element, sol);

            // update the volvars of the element
            volumeVariables_[eIdx].resize(fvGeometry.numScv());
            for (auto&& scv : scvs(fvGeometry))
            {
                volumeVariables_[eIdx][scv.index()].update(elemSol,
                                                           problem,
                                                           element,
                                                           scv);
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementVolumeVariables localView(const BoxGlobalVolumeVariables& global)
    { return ElementVolumeVariables(global); }

    const VolumeVariables& volVars(const IndexType eIdx, const IndexType scvIdx) const
    { return volumeVariables_[eIdx][scvIdx]; }

    VolumeVariables& volVars(const IndexType eIdx, const IndexType scvIdx)
    { return volumeVariables_[eIdx][scvIdx]; }

  private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    IndexType eIdx_;
    std::vector<std::vector<VolumeVariables>> volumeVariables_;
};


// Specialization when the current volume variables are not stored
template<class TypeTag>
class BoxGlobalVolumeVariables<TypeTag, /*enableGlobalVolVarCache*/false>
{
    // prev vol vars have to be a friend class in order for the assignment operator to work
    friend BoxElementVolumeVariables<TypeTag, false>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:

    void update(Problem& problem, const SolutionVector& sol)
    { problemPtr_ = &problem; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementVolumeVariables localView(const BoxGlobalVolumeVariables& global)
    { return ElementVolumeVariables(global); }

private:

    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;
};

} // end namespace

#endif
