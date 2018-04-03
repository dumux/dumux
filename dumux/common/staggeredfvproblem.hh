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
 * \ingroup Common
 * \brief Base class for all problems
 */
#ifndef DUMUX_STAGGERD_FV_PROBLEM_HH
#define DUMUX_STAGGERD_FV_PROBLEM_HH

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
#include <dune/common/rangeutilities.hh>
#else
#include <dumux/common/intrange.hh>
#endif

#include <dumux/common/properties.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux {

/*!
 * \ingroup Problems
 * \ingroup Common
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
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    /*!
     * \brief Constructor
    * \param fvGridGeometry The finite volume grid geometry
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
                auto initPriVars = asImp_().initial(scv);
                asImp_().applyInititalCellCenterSolution(sol, scv, initPriVars);
            }

            // loop over faces
            for(auto&& scvf : scvfs(fvGeometry))
            {
                auto initPriVars = asImp_().initial(scvf);
                asImp_().applyInititalFaceSolution(sol, scvf, initPriVars);
            }
        }
    }


    //! Applys the initial cell center solution
    void applyInititalCellCenterSolution(SolutionVector& sol,
                                         const SubControlVolume& scv,
                                         const PrimaryVariables& initSol) const
    {
        for(auto&& i : priVarIndices_(cellCenterIdx))
            sol[cellCenterIdx][scv.dofIndex()][i] = initSol[i];
    }

    //! Applys the initial face solution
    void applyInititalFaceSolution(SolutionVector& sol,
                                   const SubControlVolumeFace& scvf,
                                   const PrimaryVariables& initSol) const
    {
        for(auto&& i : priVarIndices_(faceIdx))
            sol[faceIdx][scvf.dofIndex()][i] = initSol[i];
    }

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for cell center dofs.
    static auto priVarIndices_(typename GET_PROP(TypeTag, DofTypeIndices)::CellCenterIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Dune::range(0, numEqCellCenter);
#else
        return IntRange(0, numEqCellCenter);
#endif
    }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for face dofs.
    static auto priVarIndices_(typename GET_PROP(TypeTag, DofTypeIndices)::FaceIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        constexpr auto numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Dune::range(numEqCellCenter, numEq);
#else
        return IntRange(numEqCellCenter, numEq);
#endif
    }

};

} // end namespace Dumux

#endif
