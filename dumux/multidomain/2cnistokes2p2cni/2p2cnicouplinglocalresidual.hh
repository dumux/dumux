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
 * \brief Extending the TwoPTwoCNILocalResidual by the required functions for
 *        a coupled application.
 */
#ifndef DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH
#define DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH

#include <dumux/porousmediumflow/nonisothermal/implicit/localresidual.hh>
#include <dumux/porousmediumflow/2p2c/implicit/indices.hh>
#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitLocalResidual
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 * \brief Extending the TwoPTwoCNILocalResidual by the required functions for
 *        a coupled application.
 */
template<class TypeTag>
class TwoPTwoCNICouplingLocalResidual : public NILocalResidual<TypeTag>
{
    typedef NILocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { useMoles = GET_PROP_VALUE(TypeTag, UseMoles) };
    enum {
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum {
        massBalanceIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx),
        contiWEqIdx = Indices::contiWEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };
    enum {
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };
    enum { phaseIdx = nPhaseIdx }; // index of the phase for the phase flux calculation
    enum { compIdx = wCompIdx}; // index of the component for the phase flux calculation

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,1> > ElementFluxVector;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*!
     * \brief Implementation of the boundary evaluation for the Darcy model
     *
     * Evaluate one part of the Dirichlet-like coupling conditions for a single
     * sub-control volume face; rest is done in the local coupling operator
     */
    void evalBoundary_()
    {
        ParentType::evalBoundary_();

        typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
        typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;
        const ReferenceElement &refElement = ReferenceElements::general(this->element_().geometry().type());

        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; scvIdx++)
        {
            // consider only SCVs on the boundary
            if (this->fvGeometry_().subContVol[scvIdx].inner)
                continue;

            // evaluate boundary conditions for the intersections of the current element
            for (const auto& intersection : intersections(this->gridView_(), this->element_()))
            {
                // handle only intersections on the boundary
                if (!intersection.boundary())
                    continue;

                // assemble the boundary for all vertices of the current face
                const int fIdx = intersection.indexInInside();
                const int numFaceVertices = refElement.size(fIdx, 1, dim);

                // loop over the single vertices on the current face
                for (int faceVertexIdx = 0; faceVertexIdx < numFaceVertices; ++faceVertexIdx)
                {
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (refElement.subEntity(fIdx, 1, faceVertexIdx, dim) != scvIdx)
                        continue;

                    const VolumeVariables &volVars = this->curVolVars_()[scvIdx];

                    // set pressure as part of the momentum coupling
                    static_assert(GET_PROP_VALUE(TypeTag, Formulation) == TwoPTwoCFormulation::pnsw,
                                  "This coupling condition is only implemented for a pnsw formulation.");
                    if (this->bcTypes_(scvIdx).isCouplingDirichlet(massBalanceIdx))
                        this->residual_[scvIdx][massBalanceIdx] = volVars.pressure(nPhaseIdx);

                    // set mass/mole fraction for transported component
                    static_assert(GET_PROP_VALUE(TypeTag, NumComponents) == 2,
                                  "This coupling condition is only implemented for two components.");
                    if (this->bcTypes_(scvIdx).isCouplingDirichlet(contiWEqIdx))
                    {
                        if (useMoles)
                            this->residual_[scvIdx][contiWEqIdx] = volVars.moleFraction(nPhaseIdx, wCompIdx);
                        else
                            this->residual_[scvIdx][contiWEqIdx] = volVars.massFraction(nPhaseIdx, wCompIdx);
                    }

                    // set temperature
                    if (this->bcTypes_(scvIdx).isCouplingDirichlet(energyEqIdx))
                        this->residual_[scvIdx][energyEqIdx] = volVars.temperature();
                }
            }
        }
    }

    /*!
     * \brief Evaluates the time derivative of the storage term
     *
     * \param storage The vector in which the result is written
     * \param scvIdx The sub-control-volume index
     */
    void evalStorageDerivative(PrimaryVariables &storage, const int scvIdx) const
    {
        PrimaryVariables result;
        this->computeStorage(result, scvIdx, false);
        Valgrind::CheckDefined(result);
        storage = result;
        this->computeStorage(result, scvIdx, true);
        Valgrind::CheckDefined(result);
        storage -= result;

        storage *= this->fvGeometry_().subContVol[scvIdx].volume
                   / this->problem_().timeManager().timeStepSize()
                   * this->curVolVars_(scvIdx).extrusionFactor();
    }
};

} // namespace Dumux

#endif // DUMUX_2P2CNI_COUPLING_LOCAL_RESIDUAL_HH
