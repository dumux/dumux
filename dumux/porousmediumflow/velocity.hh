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
 * \ingroup PorousmediumflowModels
 * \brief Velocity computation for implicit (porous media) models.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_VELOCITY_HH
#define DUMUX_POROUSMEDIUMFLOW_VELOCITY_HH

#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/flux/traits.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper structs and functions detecting if the model is compositional

template <typename T, typename ...Ts>
using MoleFractionDetector = decltype(std::declval<T>().moleFraction(std::declval<Ts>()...));

template<class T, typename ...Args>
static constexpr bool hasMoleFraction()
{ return Dune::Std::is_detected<MoleFractionDetector, T, Args...>::value; }

template <typename T, typename ...Ts>
using MassFractionDetector = decltype(std::declval<T>().massFraction(std::declval<Ts>()...));

template<class T, typename ...Args>
static constexpr bool hasMassFraction()
{ return Dune::Std::is_detected<MassFractionDetector, T, Args...>::value; }

} // end namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup PorousmediumflowModels
 * \brief Velocity computation for implicit (porous media) models.
 */
template<class GridVariables, class FluxVariables>
class PorousMediumFlowVelocity
{
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using Scalar = typename GridVariables::Scalar;
    using FluxTraits = typename Dumux::FluxTraits<FluxVariables>;
    using AdvectionType = typename FluxVariables::AdvectionType;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr int dofCodim = isBox ? dim : 0;
    static constexpr bool stationaryVelocityField = FluxTraits::hasStationaryVelocityField();

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Problem = typename GridVolumeVariables::Problem;

    static constexpr bool modelIsCompositional = Detail::hasMoleFraction<typename GridVolumeVariables::VolumeVariables, int, int>() ||
                                                 Detail::hasMassFraction<typename GridVolumeVariables::VolumeVariables, int, int>();

    static_assert(VolumeVariables::numFluidPhases() >= 1, "Velocity output only makes sense for models with fluid phases.");

public:
    static constexpr int numFluidPhases = VolumeVariables::numFluidPhases();
    using VelocityVector = std::vector<Dune::FieldVector<Scalar, dimWorld>>;

    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param gridVariables The grid variables
     */
    PorousMediumFlowVelocity(const GridVariables& gridVariables)
    : problem_(gridVariables.curGridVolVars().problem())
    , gridGeometry_(gridVariables.gridGeometry())
    , gridVariables_(gridVariables)
    {
        // set the number of scvs the vertices are connected to
        if constexpr (isBox && dim > 1)
        {
            // resize to the number of vertices of the grid
            cellNum_.assign(gridGeometry_.gridView().size(dim), 0);

            for (const auto& element : elements(gridGeometry_.gridView()))
                for (unsigned int vIdx = 0; vIdx < element.subEntities(dim); ++vIdx)
                    ++cellNum_[gridGeometry_.vertexMapper().subIndex(element, vIdx, dim)];
        }
    }

    //! Calculates the velocities for the scvs in the element.
    //! We assume the local containers to be bound to the complete stencil.
    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int phaseIdx) const
    {
        using Velocity = typename VelocityVector::value_type;

        const auto geometry = element.geometry();
        const Dune::GeometryType geomType = geometry.type();

        // the upwind term to be used for the volume flux evaluation
        auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.mobility(phaseIdx); };

        if constexpr (isBox && dim == 1)
        {
            Velocity tmpVelocity(0.0);
            tmpVelocity = (geometry.corner(1) - geometry.corner(0));
            tmpVelocity /= tmpVelocity.two_norm();

            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    continue;

                // insantiate the flux variables
                FluxVariables fluxVars;
                fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                // get the volume flux divided by the area of the
                // subcontrolvolume face in the reference element
                Scalar localArea = scvfReferenceArea_(geomType, scvf.index());
                Scalar flux = fluxVars.advectiveFlux(phaseIdx, upwindTerm) / localArea;
                const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                flux /= insideVolVars.extrusionFactor();
                tmpVelocity *= flux;

                const int eIdxGlobal = gridGeometry_.elementMapper().index(element);
                velocity[eIdxGlobal] = tmpVelocity;
            }
            return;
        }

        // get the transposed Jacobian of the element mapping
        const auto refElement = Dune::referenceElement(geometry);
        const auto& localPos = refElement.position(0, 0);
        const auto jacobianT2 = geometry.jacobianTransposed(localPos);

        if constexpr (isBox)
        {
            using ScvVelocities = Dune::BlockVector<Velocity>;
            ScvVelocities scvVelocities(fvGeometry.numScv());
            scvVelocities = 0;

            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    continue;

                // local position of integration point
                const auto localPosIP = geometry.local(scvf.ipGlobal());

                // Transformation of the global normal vector to normal vector in the reference element
                const auto jacobianT1 = geometry.jacobianTransposed(localPosIP);
                const auto globalNormal = scvf.unitOuterNormal();
                GlobalPosition localNormal(0);
                jacobianT1.mv(globalNormal, localNormal);
                localNormal /= localNormal.two_norm();

                // instantiate the flux variables
                FluxVariables fluxVars;
                fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                // get the volume flux divided by the area of the
                // subcontrolvolume face in the reference element
                Scalar localArea = scvfReferenceArea_(geomType, scvf.index());
                Scalar flux = fluxVars.advectiveFlux(phaseIdx, upwindTerm) / localArea;
                const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                flux /= insideVolVars.extrusionFactor();

                // transform the volume flux into a velocity vector
                Velocity tmpVelocity = localNormal;
                tmpVelocity *= flux;

                scvVelocities[scvf.insideScvIdx()] += tmpVelocity;
                scvVelocities[scvf.outsideScvIdx()] += tmpVelocity;
            }

            // transform vertex velocities from local to global coordinates
            for (auto&& scv : scvs(fvGeometry))
            {
                int vIdxGlobal = scv.dofIndex();

                // calculate the subcontrolvolume velocity by the Piola transformation
                Velocity scvVelocity(0);

                jacobianT2.mtv(scvVelocities[scv.indexInElement()], scvVelocity);
                scvVelocity /= geometry.integrationElement(localPos)*cellNum_[vIdxGlobal];
                // add up the wetting phase subcontrolvolume velocities for each vertex
                velocity[vIdxGlobal] += scvVelocity;
            }
        }
        else
        {
            // For the number of scvfs per facet (mpfa) we simply obtain the number of
            // corners of the first facet as prisms/pyramids are not supported here anyway
            // -> for prisms/pyramids the number of scvfs would differ from facet to facet
            static constexpr bool isMpfa = GridGeometry::discMethod == DiscretizationMethod::ccmpfa;
            const int numScvfsPerFace = isMpfa ? element.template subEntity<1>(0).geometry().corners() : 1;

            if (fvGeometry.numScvf() != element.subEntities(1)*numScvfsPerFace)
                DUNE_THROW(Dune::NotImplemented, "Velocity output for non-conforming grids");

            if (!geomType.isCube() && !geomType.isSimplex())
                DUNE_THROW(Dune::NotImplemented, "Velocity output for other geometry types than cube and simplex");

            // first we extract the corner indices for each scv for the CIV method
            // for network grids there might be multiple intersection with the same geometryInInside
            // we identify those by the indexInInside for now (assumes conforming grids at branching facets)
            // here we keep track of them
            std::vector<bool> handledScvf;
            if (dim < dimWorld)
                handledScvf.resize(element.subEntities(1), false);

            // find the local face indices of the scvfs (for conforming meshes)
            std::vector<unsigned int> scvfIndexInInside(fvGeometry.numScvf());
            int localScvfIdx = 0;
            for (const auto& intersection : intersections(gridGeometry_.gridView(), element))
            {
                if (dim < dimWorld)
                    if (handledScvf[intersection.indexInInside()])
                        continue;

                if (intersection.neighbor() || intersection.boundary())
                {
                    for (int i = 0; i < numScvfsPerFace; ++i)
                        scvfIndexInInside[localScvfIdx++] = intersection.indexInInside();

                    // for surface and network grids mark that we handled this face
                    if (dim < dimWorld)
                        handledScvf[intersection.indexInInside()] = true;
                }
            }

            std::vector<Scalar> scvfFluxes(element.subEntities(1), 0.0);
            localScvfIdx = 0;
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                    scvfFluxes[scvfIndexInInside[localScvfIdx]] += fluxVars.advectiveFlux(phaseIdx, upwindTerm)/insideVolVars.extrusionFactor();
                }
                else
                {
                    auto bcTypes = problem_.boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                        scvfFluxes[scvfIndexInInside[localScvfIdx]] += fluxVars.advectiveFlux(phaseIdx, upwindTerm)/insideVolVars.extrusionFactor();
                    }
                }

                // increment scvf counter
                localScvfIdx++;
            }

            // Correct boundary fluxes in case of Neumann conditions.
            // In this general setting, it would be very difficult to
            // calculate correct phase, i.e., volume, fluxes from arbitrary
            // Neumann conditions. We approximate the Neumann flux by the
            // flux on the opposite face. For extremely distorted grids this can
            // lead to unexpected results (but then TPFA also leads to unexpected results).
            localScvfIdx = 0;
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    auto bcTypes = problem_.boundaryTypes(element, scvf);
                    if (bcTypes.hasNeumann())
                    {
                        // for stationary velocity fields we can easily compute the correct velocity
                        // this special treatment makes sure that the velocity field is also correct on Neumann boundaries
                        // of tracer problems where the velocity field is given.
                        // (For Dirichlet boundaries no special treatment is necessary.)
                        if (stationaryVelocityField)
                        {
                            const auto flux = AdvectionType::flux(problem_, element, fvGeometry, elemVolVars,
                                                                  scvf, phaseIdx, elemFluxVarsCache);

                            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                            scvfFluxes[scvfIndexInInside[localScvfIdx]] += flux / insideVolVars.extrusionFactor();
                        }
                        else
                        {
                            // check if we have Neumann no flow, we can just use 0
                            const auto neumannFlux = problem_.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                            using NumEqVector = std::decay_t<decltype(neumannFlux)>;
                            if (Dune::FloatCmp::eq<NumEqVector, Dune::FloatCmp::CmpStyle::absolute>(neumannFlux, NumEqVector(0.0), 1e-30))
                                scvfFluxes[scvfIndexInInside[localScvfIdx]] = 0;

                            // otherwise, we try some reconstruction
                            // for cubes
                            else if (dim == 1 || geomType.isCube())
                            {
                                const auto fIdx = scvfIndexInInside[localScvfIdx];

                                if constexpr (!modelIsCompositional)
                                {
                                    // We assume that the density at the face equals the one at the cell center and reconstruct a volume flux from the Neumann mass flux.
                                    const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                                    const auto eqIdx = VolumeVariables::Indices::conti0EqIdx + phaseIdx;
                                    scvfFluxes[fIdx] += neumannFlux[eqIdx] / insideVolVars.density(phaseIdx) * scvf.area() * insideVolVars.extrusionFactor();
                                }
                                else
                                {
                                    // For compositional models, we generally can't reconstruct the volume flow from the Neumann flux (which is a species flux rather
                                    // than a phase flux here). Instead, we use the velocity of the opposing face.
                                    const auto fIdxOpposite = fIdx%2 ? fIdx-1 : fIdx+1;
                                    scvfFluxes[fIdx] = -scvfFluxes[fIdxOpposite];
                                }
                            }

                            // for simplices
                            else if (geomType.isSimplex())
                                scvfFluxes[scvfIndexInInside[localScvfIdx]] = 0;

                            else
                                DUNE_THROW(Dune::NotImplemented, "Velocity computation at Neumann boundaries for cell-centered and prism/pyramid");
                        }
                    }
                }

                // increment scvf counter
                localScvfIdx++;
            }


            Velocity refVelocity;
            // cubes: On the reference element simply average over opposite fluxes
            // note that this is equal to a corner velocity interpolation method
            if (dim == 1 || geomType.isCube())
            {
                for (int i = 0; i < dim; i++)
                    refVelocity[i] = 0.5 * (scvfFluxes[2*i + 1] - scvfFluxes[2*i]);
            }
            // simplices: Raviart-Thomas-0 interpolation evaluated at the cell center
            else if (geomType.isSimplex())
            {
                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                {
                    refVelocity[dimIdx] = -scvfFluxes[dim - 1 - dimIdx];
                    for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                        refVelocity[dimIdx] += scvfFluxes[fIdx]/(dim + 1);
                }
            }
            // 3D prism and pyramids
            else
                DUNE_THROW(Dune::NotImplemented, "Velocity computation for cell-centered and prism/pyramid");

            Velocity scvVelocity(0);
            jacobianT2.mtv(refVelocity, scvVelocity);

            scvVelocity /= geometry.integrationElement(localPos);

            int eIdxGlobal = gridGeometry_.elementMapper().index(element);

            velocity[eIdxGlobal] = scvVelocity;

        } // cell-centered
    }

private:
    // The area of a subcontrolvolume face in a reference element.
    // The 3d non-cube values have been calculated with quadrilateralArea3D
    // of boxfvelementgeometry.hh.
    static Scalar scvfReferenceArea_(Dune::GeometryType geomType, int fIdx)
    {
        if (dim == 1 || geomType.isCube())
        {
            return 1.0/(1 << (dim-1));
        }
        else if (geomType.isTriangle())
        {
            static const Scalar faceToArea[] = {0.372677996249965,
                                                0.372677996249965,
                                                0.235702260395516};
            return faceToArea[fIdx];
        }
        else if (geomType.isTetrahedron())
        {
            static const Scalar faceToArea[] = {0.102062072615966,
                                                0.102062072615966,
                                                0.058925565098879,
                                                0.102062072615966,
                                                0.058925565098879,
                                                0.058925565098879};
            return faceToArea[fIdx];
        }
        else if (geomType.isPyramid())
        {
            static const Scalar faceToArea[] = {0.130437298687488,
                                                0.130437298687488,
                                                0.130437298687488,
                                                0.130437298687488,
                                                0.150923085635624,
                                                0.1092906420717,
                                                0.1092906420717,
                                                0.0781735959970571};
            return faceToArea[fIdx];
        }
        else if (geomType.isPrism())
        {
            static const Scalar faceToArea[] = {0.166666666666667,
                                                0.166666666666667,
                                                0.166666666666667,
                                                0.186338998124982,
                                                0.186338998124982,
                                                0.117851130197758,
                                                0.186338998124982,
                                                0.186338998124982,
                                                0.117851130197758};
            return faceToArea[fIdx];
        }
        else {
            DUNE_THROW(Dune::NotImplemented, "scvf area for unknown GeometryType");
        }
    }

private:
    const Problem& problem_;
    const GridGeometry& gridGeometry_;
    const GridVariables& gridVariables_;

    std::vector<int> cellNum_;
};

} // end namespace Dumux

#endif
