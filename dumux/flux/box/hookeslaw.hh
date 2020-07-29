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
 * \ingroup BoxFlux
 * \brief Specialization of Hooke's law for the box scheme. This computes
 *        the stress tensor and surface forces resulting from mechanical deformation.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_HOOKES_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_HOOKES_LAW_HH

#include <dune/common/fmatrix.hh>

#include <dumux/flux/hookeslaw.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup BoxFlux
 * \brief Hooke's law for box scheme
 * \tparam ScalarType the scalar type for scalar physical quantities
 * \tparam GridGeometry the grid geometry
 */
template<class ScalarType, class GridGeometry>
class HookesLaw<ScalarType, GridGeometry, DiscretizationMethod::box>
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == dimWorld, "Hookes Law not implemented for network/surface grids");

public:
    //! export the type used for scalar values
    using Scalar = ScalarType;
    //! export the type used for the stress tensor
    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    //! export the type used for force vectors
    using ForceVector = typename StressTensor::row_type;
    //! state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::box;

    /*!
     * \brief Returns the force (in Newton) acting on a sub-control volume face.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static ForceVector force(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf,
                             const ElementFluxVarsCache& elemFluxVarCache)
    {
        const auto sigma = stressTensor(problem, element, fvGeometry, elemVolVars, elemFluxVarCache[scvf]);

        ForceVector scvfForce(0.0);
        sigma.mv(scvf.unitOuterNormal(), scvfForce);
        scvfForce *= Extrusion::area(scvf);

        return scvfForce;
    }

    //! assembles the stress tensor at a given integration point
    template<class Problem, class ElementVolumeVariables, class FluxVarsCache>
    static StressTensor stressTensor(const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const FluxVarsCache& fluxVarCache)
    {
        const auto& lameParams = problem.spatialParams().lameParams(element, fvGeometry, elemVolVars, fluxVarCache);

        // evaluate displacement gradient
        StressTensor gradU(0.0);
        for (int dir = 0; dir < dim; ++dir)
            for (const auto& scv : scvs(fvGeometry))
                gradU[dir].axpy(elemVolVars[scv.indexInElement()].displacement(dir), fluxVarCache.gradN(scv.indexInElement()));

        // evaluate strain tensor
        StressTensor epsilon;
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dimWorld; ++j)
                epsilon[i][j] = 0.5*(gradU[i][j] + gradU[j][i]);

        // calculate sigma
        StressTensor sigma(0.0);
        const auto traceEpsilon = trace(epsilon);
        for (int i = 0; i < dim; ++i)
        {
            sigma[i][i] = lameParams.lambda()*traceEpsilon;
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] += 2.0*lameParams.mu()*epsilon[i][j];
        }

        return sigma;
    }
};

} // end namespace Dumux

#endif
