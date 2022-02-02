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
 * \ingroup FreeflowTurbulenceModels
 * \brief Base class for all turbulence problems.
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_TURBULENCE_PROBLEM_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_TURBULENCE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/walldistance.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/turbulence/boundarytypes.hh>

namespace Dumux {
/*!
 * \ingroup FreeflowTurbulenceModels
 * \brief Base class for all turbulence problems.
 */
//! The implementation is specialized for the different discretizations
template<class TypeTag, class DiscretizationMethod> struct RANSParentProblemImpl;

//! The actual NavierStokesParentProblem
template<class TypeTag>
using RANSParentProblem =
    typename RANSParentProblemImpl<TypeTag,
                                   typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>::type;

template<class TypeTag, class DiscretizationMethod>
class RANSProblemImpl;

template<class TypeTag>
class RANSProblemImpl<TypeTag, DiscretizationMethods::FCStaggered>
: public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
public:
    using BoundaryTypes = Dumux::RANSFCBoundaryTypes<ModelTraits, ModelTraits::dim()>;

    RANSProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, couplingManager)
    { }
};

template<class TypeTag>
class RANSProblemImpl<TypeTag, DiscretizationMethods::CCTpfa>
: public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

    static constexpr bool isCoupled_ = !std::is_empty_v<CouplingManager>;
public:
    //! Export spatial parameter type
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    //! Export the boundary types.
    using BoundaryTypes = Dumux::RANSCCBoundaryTypes<ModelTraits, ModelTraits::numEq()>;
    /*!
     * \brief Constructor, passing the spatial parameters.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param spatialParams The spatial parameter class
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    RANSProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    std::shared_ptr<SpatialParams> spatialParams,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, couplingManager)
    , spatialParams_(spatialParams)
    {
        isFlatWallBounded_ = getParam<bool>("RANS.IsFlatWallBounded", false);
        twoEqTurbulenceModelName_ = getParam<std::string>("RANS.TwoEqTurbulenceModelName");
        checkForWalls_();
        manageWallInformation_();
    }

    /*!
     * \brief Constructor, constructing the spatial parameters.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    RANSProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : RANSProblemImpl(gridGeometry,
                      couplingManager,
                      std::make_shared<SpatialParams>(gridGeometry),
                      paramGroup)
    { }

    /*!
     * \brief Returns the normal velocity at a given sub control volume face.
     */
    Scalar stressTensorScalarProduct(const Element& element,
                                     const SubControlVolume& scv) const
    {
        if constexpr (isCoupled_)
            return this->couplingManager().stressTensorScalarProduct(element, scv);
        else
            return asImp_().stressTensorScalarProductAtPos(scv.ipGlobal());
    }

    /*!
     * \brief Returns the velocity at a given position.
     */
    Scalar stressTensorScalarProductAtPos(const GlobalPosition&) const
    { DUNE_THROW(Dune::NotImplemented, "stressTensorScalarProductAtPos not implemented"); }

    /*!
     * \brief Returns the normal velocity at a given sub control volume face.
     */
    Scalar vorticityTensorScalarProduct(const Element& element,
                                        const SubControlVolume& scv) const
    {
        if constexpr (isCoupled_)
            return this->couplingManager().vorticityTensorScalarProduct(element, scv);
        else
            return asImp_().vorticityTensorScalarProductAtPos(scv.ipGlobal());
    }

    /*!
     * \brief Returns the velocity at a given position.
     */
    Scalar vorticityTensorScalarProductAtPos(const GlobalPosition&) const
    { DUNE_THROW(Dune::NotImplemented, "vorticityTensorScalarProductAtPos not implemented");    }

    /*!
     * \brief Returns the name of the turbulence model
     */
    std::string twoEqTurbulenceModelName() const
    { return twoEqTurbulenceModelName_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return *spatialParams_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

private:

    void checkForWalls_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (asImp_().boundaryTypes(element, scvf).hasWall())
                    return;
            }
        }
        // If reached, no walls were found, throw exception.
        DUNE_THROW(Dune::InvalidStateException, "No walls are are specified with the setWall() function");
    }

    /*!
     * \brief Use the boundary search algorithm to find the shortest distance to a wall for each element
     *
     *  Also store the wall element's index, and its direction in the case of flat wall bounded problems
     */
    void manageWallInformation_()
    {
        WallDistance wallInformation(this->gridGeometry(), WallDistance<GridGeometry>::atElementCenters,
            [this] (const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
            { return asImp_().boundaryTypes(fvGeometry.element(), scvf).hasWall(); });

        spatialParams().setWallDistance(wallInformation.wallDistance());
        if (isFlatWallBounded_)
            spatialParams().setWallData(wallInformation.wallData(), this->gridGeometry());
    }


protected:
    // material properties of the porous medium
    std::shared_ptr<SpatialParams> spatialParams_;
    bool isFlatWallBounded_;
    std::string twoEqTurbulenceModelName_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

/*!
 * \ingroup RANSModels
 * \brief RANS problem class
 *
 * Inherit from this problem to implement RANS problems
 */
template<class TypeTag>
using RANSProblem = RANSProblemImpl<
            TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>;
} // end namespace Dumux

#endif
