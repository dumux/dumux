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
 * \brief Base class for all problems
 */
#ifndef DUMUX_STAGGERD_FV_PROBLEM_HH
#define DUMUX_STAGGERD_FV_PROBLEM_HH

#include <dumux/discretization/staggered/properties.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux
{
/*!
 * \ingroup Problems
 * \brief Base class for all finite-volume problems
 *
 * \note All quantities (regarding the units) are specified assuming a
 *       three-dimensional world. Problems discretized using 2D grids
 *       are assumed to be extruded by \f$1 m\f$ and 1D grids are assumed
 *       to have a cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class StaggeredFVProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using PointSourceHelper = typename GET_PROP_TYPE(TypeTag, PointSourceHelper);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;


    // using GridAdaptModel = ImplicitGridAdapt<TypeTag, adaptiveGrid>;

public:
    /*!
     * \brief Constructor
     *
     * \param gridView The simulation's idea about physical space
     */
    StaggeredFVProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    { }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "an initial() or an initialAtPos() method.");
    }

    /*!
     * \brief Evaluate the initial value for
     * an element (for cell-centered models)
     * or vertex (for box / vertex-centered models)
     *
     * \param entity The dof entity (element or vertex)
     */
    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const
    {
        // static_assert(int(Entity::codimension) == 0 || int(Entity::codimension) == dim, "Entity must be element or vertex");
        return asImp_().initialAtPos(entity.center());
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     *
    */
    void applyInitialSolution(SolutionVector& sol) const
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            // loop over sub control volumes
            for (auto&& scv : scvs(fvGeometry))
            {
                // let the problem do the dirty work of nailing down
                // the initial solution.
                auto initPriVars = asImp_().initial(scv)[cellCenterIdx];
                auto dofIdxGlobal = scv.dofIndex();
                sol[cellCenterIdx][dofIdxGlobal] += initPriVars;
            }

            // loop over faces
            for(auto&& scvf : scvfs(fvGeometry))
            {
                auto initPriVars = asImp_().initial(scvf)[faceIdx][scvf.directionIndex()];
                sol[faceIdx][scvf.dofIndex()] = initPriVars;
            }
        }
    }

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

};

} // end namespace Dumux

#endif
