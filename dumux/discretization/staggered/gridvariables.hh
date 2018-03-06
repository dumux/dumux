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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVariables
 */
#ifndef DUMUX_STAGGERED_GRID_VARIABLES_HH
#define DUMUX_STAGGERED_GRID_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/fvgridvariables.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for cell center of face specific auxiliary GridVariables classes.
 *        Provides a common interface and a pointer to the actual grid variables.
 */
template<class ActualGridVariables>
class StaggeredGridVariablesView
{
    using GridVolumeVariables = typename ActualGridVariables::GridVolumeVariables;
    using GridFaceVariables = typename ActualGridVariables::GridFaceVariables;
    using GridFluxVariablesCache = typename ActualGridVariables::GridFluxVariablesCache;
    using FVGridGeometry = typename ActualGridVariables::FVGridGeometry;

public:
    explicit StaggeredGridVariablesView(ActualGridVariables* gridVariables)
    : gridVariables_(gridVariables) {}

    //! return the flux variables cache
    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridVariables_->gridFluxVarsCache(); }

    //! return the flux variables cache
    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridVariables_->gridFluxVarsCache(); }

    //! return the current volume variables
    const GridVolumeVariables& curGridVolVars() const
    { return gridVariables_->curGridVolVars(); }

    //! return the current volume variables
    GridVolumeVariables& curGridVolVars()
    { return gridVariables_->curGridVolVars(); }

    //! return the volume variables of the previous time step (for instationary problems)
    const GridVolumeVariables& prevGridVolVars() const
    { return gridVariables_->prevGridVolVars(); }

    //! return the volume variables of the previous time step (for instationary problems)
    GridVolumeVariables& prevGridVolVars()
    { return gridVariables_->prevGridVolVars(); }

    //! return the current face variables
    const GridFaceVariables& curGridFaceVars() const
    { return gridVariables_->curGridFaceVars(); }

    //! return the previous face variables
    const GridFaceVariables& prevGridFaceVars() const
    { return gridVariables_->prevGridFaceVars(); }

    //! return the current face variables
    GridFaceVariables& curGridFaceVars()
    { return gridVariables_->curGridFaceVars(); }

    //! return the previous face variables
    GridFaceVariables& prevGridFaceVars()
    { return gridVariables_->prevGridFaceVars(); }

    //! return the fv grid geometry
    const FVGridGeometry& fvGridGeometry() const
    { return (*gridVariables_->fvGridGeometry_);    }

    // return the actual grid variables
    const ActualGridVariables& gridVariables() const
    { return *gridVariables_; }

    // return the actual grid variables
    ActualGridVariables& gridVariables()
    { return *gridVariables_; }

protected:
    ActualGridVariables* gridVariables_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Cell center specific auxiliary GridVariables classes.
 *        Required for the Dumux multi-domain framework.
 */
template<class ActualGridVariables>
class CellCenterGridVariablesView : public StaggeredGridVariablesView<ActualGridVariables>
{
    using ParentType = StaggeredGridVariablesView<ActualGridVariables>;
public:
    using ParentType::ParentType;

    //! initialize all variables (stationary case)
    template<class SolVector>
    void init(const SolVector& curSol)
    {
        this->curGridVolVars().update(this->fvGridGeometry(), curSol);
        this->gridFluxVarsCache().update(this->fvGridGeometry(), this->curGridVolVars(), curSol, true);
    }

    //! initialize all variables (instationary case)
    template<class SolVector>
    void init(const SolVector& curSol, const SolVector& initSol)
    {
        this->curGridVolVars().update(this->fvGridGeometry(), curSol);
        this->gridFluxVarsCache().update(this->fvGridGeometry(), this->curGridVolVars(), curSol, true);
        this->prevGridVolVars().update(this->fvGridGeometry(), initSol);
    }

    //! update the volume variables and the flux variables cache
    template<class SolVector>
    void update(const SolVector& curSol)
    {
        this->curGridVolVars().update(this->fvGridGeometry(), curSol);
        this->gridFluxVarsCache().update(this->fvGridGeometry(), this->curGridVolVars(), curSol);
    }

    //! resets state to the one before time integration
    template<class SolVector>
    void resetTimeStep(const SolVector& sol)
    {
        this->curGridVolVars() = this->prevGridVolVars();
        this->gridFluxVarsCache().update(this->fvGridGeometry(), this->curGridVolVars(), sol);
    }
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face specific auxiliary GridVariables classes.
 *        Required for the Dumux multi-domain framework.
 */
template<class ActualGridVariables>
class FaceGridVariablesView : public StaggeredGridVariablesView<ActualGridVariables>
{
    using ParentType = StaggeredGridVariablesView<ActualGridVariables>;
public:
    using ParentType::ParentType;

    //! initialize all variables (stationary case)
    template<class SolVector>
    void init(const SolVector& curSol)
    {
        this->curGridFaceVars().update(this->fvGridGeometry(), curSol);
    }

    //! initialize all variables (instationary case)
    template<class SolVector>
    void init(const SolVector& curSol, const SolVector& initSol)
    {
        this->curGridFaceVars().update(this->fvGridGeometry(), curSol);
        this->prevGridFaceVars().update(this->fvGridGeometry(), initSol);
    }

    //! update the face variables
    template<class SolVector>
    void update(const SolVector& curSol)
    {
        this->curGridFaceVars().update(this->fvGridGeometry(), curSol);
    }

    //! resets state to the one before time integration
    template<class SolVector>
    void resetTimeStep(const SolVector& sol)
    {
        this->curGridFaceVars() = this->prevGridFaceVars();
    }
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class storing data associated to scvs and scvfs
 */
template<class TypeTag>
class StaggeredGridVariables : public FVGridVariables<TypeTag>
{
    using ParentType = FVGridVariables<TypeTag>;
    using ThisType = StaggeredGridVariables<TypeTag>;

    friend class StaggeredGridVariablesView<ThisType>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    //! Constructor
    StaggeredGridVariables(std::shared_ptr<const Problem> problem,
                           std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(problem, fvGridGeometry)
    , fvGridGeometry_(fvGridGeometry)
    , curGridFaceVariables_(*problem)
    , prevGridFaceVariables_(*problem)
    {}

    using CellCenterGridVariablesType = CellCenterGridVariablesView<ThisType>;
    using FaceGridVariablesType = FaceGridVariablesView<ThisType>;

    using FVGridGeometryTuple = std::tuple< CellCenterGridVariablesType, FaceGridVariablesType >;

    //! initialize with the complete solution vector
    void init(const SolutionVector& curSol)
    {
        ParentType::init(curSol);
        curGridFaceVariables_.update(*fvGridGeometry_, curSol);
    }

    //! initialize with the complete solution vector
    void init(const SolutionVector& curSol, const SolutionVector& initSol)
    {
        ParentType::init(curSol, initSol);
        curGridFaceVariables_.update(*fvGridGeometry_, curSol);
        prevGridFaceVariables_.update(*fvGridGeometry_, initSol);
    }

    //! Sets the current state as the previous for next time step
    //! this has to be called at the end of each time step
    void advanceTimeStep()
    {
        ParentType::advanceTimeStep();
        prevGridFaceVariables_ = curGridFaceVariables_;
    }

    //! return the current face variables
    const GridFaceVariables& curGridFaceVars() const
    { return curGridFaceVariables_; }

    //! return the previous face variables
    const GridFaceVariables& prevGridFaceVars() const
    { return prevGridFaceVariables_; }

    //! return the current face variables
    GridFaceVariables& curGridFaceVars()
    { return curGridFaceVariables_; }

    //! return the previous face variables
    GridFaceVariables& prevGridFaceVars()
    { return prevGridFaceVariables_; }

    //! Returns a pointer the cell center specific auxiliary class. Required for the multi-domain FVAssembler's ctor.
    std::unique_ptr<CellCenterGridVariablesView<ThisType>> cellCenterGridVariablesPtr()
    {
        return std::make_unique<CellCenterGridVariablesView<ThisType>>(this);
    }

    //! Returns a pointer the face specific auxiliary class. Required for the multi-domain FVAssembler's ctor.
    std::unique_ptr<FaceGridVariablesView<ThisType>> faceGridVariablesPtr()
    {
        return std::make_unique<FaceGridVariablesView<ThisType>>(this);
    }

    //! Return a copy of the cell center specific auxiliary class.
    CellCenterGridVariablesView<ThisType> cellCenterGridVariables() const
    {
        return CellCenterGridVariablesView<ThisType>(this);
    }

    //! Return a copy of the face specific auxiliary class.
    FaceGridVariablesView<ThisType> faceGridVariables() const
    {
        return FaceGridVariablesView<ThisType>(this);
    }

private:

    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    GridFaceVariables curGridFaceVariables_;
    GridFaceVariables prevGridFaceVariables_;
};

} // end namespace Dumux

#endif
