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
 * \brief Element-wise calculation the local Jacobian for the
 *        linear elastic model in the fully implicit scheme.
 */
#ifndef DUMUX_ELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_ELASTIC_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 *
 * \ingroup ElasticBoxModel
 * \ingroup ImplicitLocalResidual
 * \brief Calculate the local Jacobian for the linear
 *        elasticity model
 *
 * This class is used to fill the gaps in BoxLocalResidual for
 * the linear elasticity model.
 */
template<class TypeTag>
class ElasticLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

public:
    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        within a finite volume.
     *
     *        \param storage The storage of a quantity in the sub-control volume
     *        \param scvIdx The index of the considered face of the sub-control volume
     *        \param usePrevSol Evaluate function with solution of current or previous time step
     */
    PrimaryVariables computeStorage(const SubControlVolume& scv, const VolumeVariables& volVars) const
    {
        // quasistationary conditions assumed
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluate the stress across a face of a sub-control
     *        volume.
     *
     *        \param flux The stress over the SCV (sub-control-volume) face
     *        \param fIdx The index of the considered face of the sub control volume
     *        \param onBoundary A boolean variable to specify whether the flux variables
     *               are calculated for interior SCV faces or boundary faces, default=false
     */
    PrimaryVariables computeFlux(const SubControlVolumeFace& scvFace) const
    {
        FluxVariables fluxVars;
        fluxVars.initAndComputeFluxes(this->problem_(), this->element_(), scvFace);
        return fluxVars.stressVector();
    }

    /*!
     * \brief Calculate the source term of the equation
     *        \param source The source/sink in the SCV is the gravity term in the momentum balance
     *        \param scvIdx The index of the vertex of the sub control volume
     *
     */
    PrimaryVariables computeSource(const SubControlVolume& scv)
    {
        PrimaryVariables source(0.0);

        source += ParentType::computeSource(scv);

        // gravity term of the solid matrix in the momentum balance
        DimVector gravityTerm(0.0);
        gravityTerm = this->problem_().gravity();
        gravityTerm *= this->problem_().model().curVolVars(scv).rockDensity();

        for (int i = 0; i < dim; ++i)
          source[Indices::momentum(i)] += gravityTerm[i];

        return source;
    }

};

} // end namespace Dumux

#endif // DUMUX_ELASTIC_LOCAL_RESIDUAL_HH
