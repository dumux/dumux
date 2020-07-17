// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
#ifndef DUMUX_NAVIERSTOKES_MASS_AND_ENERGY_PROBLEM_HH
#define DUMUX_NAVIERSTOKES_MASS_AND_ENERGY_PROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/massandenergy/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem base class.
 *
 * This implements gravity (if desired) and a some other functions.
 *
 */
template<class TypeTag>
class NavierStokesMassAndEnergyProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:

    //! export the boundary types
    using BoundaryTypes = NavierStokesMassAndEnergyBoundaryTypes<ModelTraits::numEq(), typename ModelTraits::Indices>;

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesMassAndEnergyProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<CouplingManager> couplingManager,
                                     const std::string& paramGroup = "")
    : NavierStokesMassAndEnergyProblem(gridGeometry, paramGroup)
    {
         couplingManager_ = couplingManager;
    }

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesMassAndEnergyProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                     const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup) {}


    /*!
     * \brief Returns the normal velocity at a given sub control volume face.
     */
    Scalar faceVelocity(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const SubControlVolumeFace& scvf) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
        {
            auto v = asImp_().velocityAtPos(scvf.ipGlobal());
            v *= scvf.unitOuterNormal();
            return v;
        }
        else
            return couplingManager_->faceVelocity(element, scvf);
    }

    /*!
     * \brief Returns the velocity at a given position.
     */
    GlobalPosition velocityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "velocityAtPos not implemented");
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     * \note  Overload this if a fixed density shall be prescribed.
     */
    Scalar density(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const SubControlVolumeFace& scvf) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().densityAtPos(scvf.ipGlobal());
        else
            return couplingManager_->density(element, fvGeometry, scvf);
    }

    /*!
     * \brief Returns the density at a given position.
     */
    Scalar densityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "densityAtPos not implemented");
    }

private:

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::shared_ptr<CouplingManager> couplingManager_;


};

} // end namespace Dumux

#endif
