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
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePSpatialParams : public FVSpatialParamsOneP< typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                      typename GET_PROP_TYPE(TypeTag, Scalar),
                                                      OnePSpatialParams<TypeTag> >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = OnePSpatialParams<TypeTag>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    {}

    //! Return the permeability at a given position
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPoss) const
    { return permeability_; }

    //! Return the porosity for a sub-control volume
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        // let the cooupling manager evaluate divergence of displacement
        const auto divU = couplingManager().computeDivU(element, scv.center());
        return porosity_*divU;
    }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<const CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace

#endif
