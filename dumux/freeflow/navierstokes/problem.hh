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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesProblem
 */
#ifndef DUMUX_NAVIERSTOKES_PROBLEM_HH
#define DUMUX_NAVIERSTOKES_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/methods.hh>
#include "model.hh"

namespace Dumux
{


//! The implementation is specialized for the different discretizations
template<class TypeTag, DiscretizationMethods DM> struct NavierStokesParentProblemImpl;

template<class TypeTag>
struct NavierStokesParentProblemImpl<TypeTag, DiscretizationMethods::Staggered>
{
    using type = StaggeredFVProblem<TypeTag>;
};

//! The actual NavierStokesParentProblem
template<class TypeTag>
using NavierStokesParentProblem =
      typename NavierStokesParentProblemImpl<TypeTag,
      GET_PROP_VALUE(TypeTag, DiscretizationMethod)>::type;

/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem base class.
 *
 * This implements gravity (if desired) and a function returning the temperature.
 * Includes a specialized method used only by the staggered grid discretization.
  *
 */
template<class TypeTag>
class NavierStokesProblem : public NavierStokesParentProblem<TypeTag>
{
    using ParentType = NavierStokesParentProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
      };
    // TODO: dim or dimWorld appropriate here?
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! The constructor sets the gravity, if desired by the user.
    NavierStokesProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry),
          gravity_(0)
    {
        if (getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>Problem.EnableGravity</tt> parameter is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    //! Applys the initial face solution (velocities on the faces). Specialization for staggered grid discretization.
    template <class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, DiscretizationMethod) == DiscretizationMethods::Staggered, void>::type
    applyInititalFaceSolution(SolutionVector& sol,
                              const SubControlVolumeFace& scvf,
                              const PrimaryVariables& initSol) const
    {
        typename GET_PROP(TypeTag, DofTypeIndices)::FaceIdx faceIdx;
        const auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        sol[faceIdx][scvf.dofIndex()][numEqCellCenter] = initSol[Indices::velocity(scvf.directionIndex())];
    }

private:

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;
};

}

#endif
