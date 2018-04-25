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

#include <dune/common/rangeutilities.hh>

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
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = FVGridGeometry::faceIdx();

public:
    /*!
     * \brief Constructor
     * \param fvGridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    StaggeredFVProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup)
    { }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume (-face).
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elementFaceVars All face variables for the element
     * \param e The geometrical entity on which the source shall be applied (scv or scvf)
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables, class ElementFaceVariables, class Entity>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementFaceVariables& elementFaceVars,
                       const Entity &e) const
    {
        // forward to solution independent, fully-implicit specific interface
        return asImp_().sourceAtPos(e.center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     * \param elemFaceVars All face variables for the element
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        // forward it to the interface with only the global position
        return asImp_().neumannAtPos(scvf.ipGlobal());
    }

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
     * \brief Evaluate the initial value for an element (for cell-centered primary variables)
     * or face (for velocities)
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
                asImp_().applyInitialCellCenterSolution(sol, scv, initPriVars);
            }

            // loop over faces
            for(auto&& scvf : scvfs(fvGeometry))
            {
                auto initPriVars = asImp_().initial(scvf);
                asImp_().applyInitialFaceSolution(sol, scvf, initPriVars);
            }
        }
    }


    //! Applys the initial cell center solution
    void applyInitialCellCenterSolution(SolutionVector& sol,
                                         const SubControlVolume& scv,
                                         const PrimaryVariables& initSol) const
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        constexpr auto numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
        constexpr auto offset = numEq - numEqCellCenter;

        for(auto&& i : priVarIndices_(cellCenterIdx))
            sol[cellCenterIdx][scv.dofIndex()][i] = initSol[i + offset];
    }

    //! Applys the initial face solution
    void applyInitialFaceSolution(SolutionVector& sol,
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
    static auto priVarIndices_(typename FVGridGeometry::DofTypeIndices::CellCenterIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        return Dune::range(0, numEqCellCenter);
    }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for face dofs.
    static auto priVarIndices_(typename FVGridGeometry::DofTypeIndices::FaceIdx)
    {
        constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
        return Dune::range(0, numEqFace);
    }

};

} // end namespace Dumux

#endif
