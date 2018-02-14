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
 * \ingroup SweModel
 * \copydoc Dumux::SweProblem
 */
#ifndef DUMUX_SWE_PROBLEM_HH
#define DUMUX_SWE_PROBLEM_HH

#include <dumux/common/fvproblem.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include "model.hh"

namespace Dumux
{


/*!
 * \ingroup SweModel
 * \brief Swe problem base class.
 *
 *
 */
template<class TypeTag>
class SweProblem : public FVProblem<TypeTag>
{
    using ParentType = FVPProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Element = typename GridView::template Codim<0>::Entity;

    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
      };
    // TODO: dim or dimWorld appropriate here?
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    Scalar getH(const element)
    {
        return asImp_().getH();
    }

    Scalar getU(const element)
    {
        return asImp_().getU();
    }

    Scalar getVH(const element)
    {
        return asImp_().getV();
    }

    Scalar getVH(const element)
    {
        return asImp_().getZ();
    }

    Scalar getFrictionH(const element)
    {
        return asImp_().getFrictionH();
    }

    Scalar getFrictionUstarH(const element)
    {
        return asImp_().getFrictionH();
    }

    Scalar getFrictionUstarH(const element)
    {
        return asImp_().getGravity();
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
