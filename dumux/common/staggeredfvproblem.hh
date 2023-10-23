// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \ingroup StaggeredDiscretization
 * \brief Base class for all staggered fv problems
 */
#ifndef DUMUX_STAGGERD_FV_PROBLEM_HH
#define DUMUX_STAGGERD_FV_PROBLEM_HH

#include <dune/common/rangeutilities.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

/*!
 * \ingroup Core
 * \ingroup StaggeredDiscretization
 * \brief Base class for all staggered finite-volume problems
 *
 * \note All quantities (regarding the units) are specified assuming a
 *       three-dimensional world. Problems discretized using 2D grids
 *       are assumed to be extruded by \f$1 m\f$ and 1D grids are assumed
 *       to have a cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class StaggeredFVProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto cellCenterIdx = GridGeometry::cellCenterIdx();
    static constexpr auto faceIdx = GridGeometry::faceIdx();

    static constexpr auto numEqCellCenter = getPropValue<TypeTag, Properties::NumEqCellCenter>();
    static constexpr auto numEqFace = getPropValue<TypeTag, Properties::NumEqFace>();

public:
    /*!
     * \brief Constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    StaggeredFVProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                       const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    { }

    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     * \param pvIdx The primary variable index
     */
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    { return false; }

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
        return this->asImp_().sourceAtPos(e.center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
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
        return this->asImp_().neumannAtPos(scvf.ipGlobal());
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
        return this->asImp_().initialAtPos(entity.center());
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     *
    */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol[cellCenterIdx].resize(this->gridGeometry().numCellCenterDofs());
        sol[faceIdx].resize(this->gridGeometry().numFaceDofs());

        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);

            // loop over sub control volumes
            for (auto&& scv : scvs(fvGeometry))
            {
                // let the problem do the dirty work of nailing down
                // the initial solution.
                auto initPriVars = this->asImp_().initial(scv);
                this->asImp_().applyInitialCellCenterSolution(sol, scv, initPriVars);
            }

            // loop over faces
            for(auto&& scvf : scvfs(fvGeometry))
            {
                auto initPriVars = this->asImp_().initial(scvf);
                this->asImp_().applyInitialFaceSolution(sol, scvf, initPriVars);
            }
        }
    }


    //! Applies the initial cell center solution
    template<class SolutionVector>
    void applyInitialCellCenterSolution(SolutionVector& sol,
                                        const SubControlVolume& scv,
                                        const PrimaryVariables& initSol) const
    {
        // while the container within the actual solution vector holds numEqCellCenter
        // elements, we need to specify an offset to get the correct entry of the initial solution
        static constexpr auto offset = PrimaryVariables::dimension - numEqCellCenter;

        for(int pvIdx = 0; pvIdx < numEqCellCenter; ++pvIdx)
            sol[cellCenterIdx][scv.dofIndex()][pvIdx] = initSol[pvIdx + offset];
    }

    //! Applies the initial face solution
    template<class SolutionVector>
    void applyInitialFaceSolution(SolutionVector& sol,
                                  const SubControlVolumeFace& scvf,
                                  const PrimaryVariables& initSol) const
    {
        for(int pvIdx = 0; pvIdx < numEqFace; ++pvIdx)
            sol[faceIdx][scvf.dofIndex()][pvIdx] = initSol[pvIdx];
    }

};

} // end namespace Dumux

#endif
