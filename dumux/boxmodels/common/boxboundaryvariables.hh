/*****************************************************************************
 *   Copyright (C) 2011 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the data which is required to calculate
 *        the component fluxes of the fluid phase over the boundary of a finite volume.
 */
#ifndef DUMUX_BOX_BOUNDARY_VARIABLES_HH
#define DUMUX_BOX_BOUNDARY_VARIABLES_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \brief This template class contains data which is required to
 *        calculate the component fluxes of the fluid phase over the boundary of a
 *        finite volume a box model.
 *
 * The class is derived from the actual FluxVariables class of the model and
 * uses the integration point of the boundary face for the evaluation of the quantities.
 * This means e.g. concentration gradients, diffusion coefficient, mass fraction at
 * the integration point.
 */
template <class TypeTag>
class BoxBoundaryVariables : public GET_PROP_TYPE(TypeTag, PTAG(FluxVariables))
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename FVElementGeometry::BoundaryFace BoundaryFace;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };

public:
    BoxBoundaryVariables(const Problem &problem,
                         const Element &element,
                         const FVElementGeometry &elemGeom,
                         int boundaryFaceIdx,
                         const ElementVolumeVariables &elemDat)
        : ParentType(problem, element, elemGeom, boundaryFaceIdx, elemDat, /*onBoundary=*/ true),
          fvGeom_(elemGeom)
    {
        // evaluate variables at the integration point of the boundary face
        boundaryFaceIdx_ = boundaryFaceIdx;
        const BoundaryFace &boundaryFace = elemGeom.boundaryFace[boundaryFaceIdx_];

        ParentType::calculateValues_(problem, element, boundaryFace, elemDat,  /*onBoundary=*/ true);
    };

    // CAREFUL: use this method to address boundary faces
    // DO NOT USE face() from the parent flux variables!!!
    const BoundaryFace &boundaryFace() const
    { return fvGeom_.boundaryFace[boundaryFaceIdx_]; };

private:
    FVElementGeometry fvGeom_;
    int boundaryFaceIdx_;
};

} // end namespace

#endif
