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
 *
 * \brief This file contains the calculation of all the fluxes over the surface of the
 * finite volume that make up the volume, the mass and the momentum balance
 * for the two-phase linear-elastic model.
 *
 * This means pressure, concentration and solid-displacement gradients, phase densities at
 * the integration point, etc.
 *
 * This class inherits from the two-phase model FluxVariables
 */
#ifndef DUMUX_DECOUPLED_ELASTIC_FLUX_VARIABLES_HH
#define DUMUX_DECOUPLED_ELASTIC_FLUX_VARIABLES_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include "dumux/geomechanics/el2p/properties.hh"

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(SpatialParams);
}
/*!
 * \ingroup ElTwoPBoxModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the fluxes over the surface of the
 *           finite volume that make up the volume, the mass and the momentum balance
 *           for the two-phase linear-elastic model.
 *
 * This means pressure, concentration and solid-displacement gradients, phase densities at
 * the integration point, etc.
 *
 */
template<class TypeTag>
class DecoupledElasticFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
                enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    enum {numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:
    /*
     * \brief The old constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    DUNE_DEPRECATED_MSG("FluxVariables now have to be default constructed and updated.")
    DecoupledElasticFluxVariables(const Problem &problem,
                        const Element &element,
                        const FVElementGeometry &fvGeometry,
                        int fIdx,
                        const ElementVolumeVariables &elemVolVars,
                        const bool onBoundary = false)
    {}

    /*!
     * \brief Default constructor
     * \note This can be removed when the deprecated constructor is removed.
     */
    DecoupledElasticFluxVariables() = default;

    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
    }

public:

    const SCVFace &face() const
    {
        return this->fvGeometry_().subContVolFace[this->faceIdx_];
    }

protected:

};

} // end namespace

#endif
