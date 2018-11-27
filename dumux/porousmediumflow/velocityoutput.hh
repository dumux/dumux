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
 * \brief Velocity output for porous media models
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH
#define DUMUX_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH

#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/velocityoutput.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \brief Velocity output policy for implicit (porous media) models
 */
template<class GridVariables, class FluxVariables>
class PorousMediumFlowVelocityOutput : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using Scalar = typename GridVariables::Scalar;

    // TODO should be possible to get this
    using Problem = typename std::decay_t<decltype(std::declval<GridVolumeVariables>().problem())>;
    using BoundaryTypes = typename std::decay_t<decltype(std::declval<GridVolumeVariables>().problem()
                                                         .boundaryTypes(std::declval<Element>(), std::declval<SubControlVolumeFace>()))>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr int dofCodim = isBox ? dim : 0;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ReferenceElements = Dune::ReferenceElements<typename GridView::ctype, dim>;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param gridVariables The grid variables
     */
    PorousMediumFlowVelocityOutput(const GridVariables& gridVariables)
    : problem_(gridVariables.curGridVolVars().problem())
    , fvGridGeometry_(gridVariables.fvGridGeometry())
    , gridVariables_(gridVariables)
    {
        // check, if velocity output can be used (works only for cubes so far)
        enableOutput_ = getParamFromGroup<bool>(problem_.paramGroup(), "Vtk.AddVelocity");
        if (enableOutput_)
        {
            // set the number of scvs the vertices are connected to
            if (isBox && dim > 1)
            {
                // resize to the number of vertices of the grid
                cellNum_.assign(fvGridGeometry_.gridView().size(dim), 0);

                for (const auto& element : elements(fvGridGeometry_.gridView()))
                    for (unsigned int vIdx = 0; vIdx < element.subEntities(dim); ++vIdx)
                        ++cellNum_[fvGridGeometry_.vertexMapper().subIndex(element, vIdx, dim)];
            }
        }
    }

    //! returns whether or not velocity output is enabled
    bool enableOutput() const override { return enableOutput_; }

    //! returns the phase name of a given phase index
    std::string phaseName(int phaseIdx) const override { return FluidSystem::phaseName(phaseIdx); }

    //! returns the number of phases
    int numPhases() const override { return VolumeVariables::numPhases(); }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    void calculateVelocity(VelocityVector& velocity,
                           const ElementVolumeVariables& elemVolVars,
                           const FVElementGeometry& fvGeometry,
                           const Element& element,
                           int phaseIdx) const override
    {
        using Velocity = typename VelocityVector::value_type;

        if (!enableOutput_) return;

        const auto geometry = element.geometry();
        const Dune::GeometryType geomType = geometry.type();

        // bind the element flux variables cache
        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        // the upwind term to be used for the volume flux evaluation
        auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.mobility(phaseIdx); };

        if(isBox && dim == 1)
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
                flux /= problem_.extrusionFactor(element,
                                                 fvGeometry.scv(scvf.insideScvIdx()),
                                                 elementSolution(element, elemVolVars, fvGeometry));
                tmpVelocity *= flux;

                const int eIdxGlobal = fvGridGeometry_.elementMapper().index(element);
                velocity[eIdxGlobal] = tmpVelocity;
            }
            return;
        }

        // get the transposed Jacobian of the element mapping
        const auto referenceElement = ReferenceElements::general(geomType);
        const auto& localPos = referenceElement.position(0, 0);
        const auto jacobianT2 = geometry.jacobianTransposed(localPos);

        if(isBox)
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
                flux /= problem_.extrusionFactor(element,
                                                 fvGeometry.scv(scvf.insideScvIdx()),
                                                 elementSolution(element, elemVolVars, fvGeometry));

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
            static constexpr bool isMpfa = FVGridGeometry::discMethod == DiscretizationMethod::ccmpfa;
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
            for (const auto& intersection : intersections(fvGridGeometry_.gridView(), element))
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
                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                    scvfFluxes[scvfIndexInInside[localScvfIdx]] += fluxVars.advectiveFlux(phaseIdx, upwindTerm)
                                                                   /problem_.extrusionFactor(element,
                                                                                             fvGeometry.scv(scvf.insideScvIdx()),
                                                                                             elementSolution(element, elemVolVars, fvGeometry));
                }
                else
                {
                    auto bcTypes = problemBoundaryTypes_(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                        scvfFluxes[scvfIndexInInside[localScvfIdx]] += fluxVars.advectiveFlux(phaseIdx, upwindTerm)
                                                                       /problem_.extrusionFactor(element,
                                                                                                 fvGeometry.scv(scvf.insideScvIdx()),
                                                                                                 elementSolution(element, elemVolVars, fvGeometry));
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
                    auto bcTypes = problemBoundaryTypes_(element, scvf);
                    if (bcTypes.hasNeumann())
                    {
                        // check if we have Neumann no flow, we can just use 0
                        const auto neumannFlux = problem_.neumann(element, fvGeometry, elemVolVars, scvf);
                        using NumEqVector = std::decay_t<decltype(neumannFlux)>;
                        if (Dune::FloatCmp::eq<NumEqVector, Dune::FloatCmp::CmpStyle::absolute>(neumannFlux, NumEqVector(0.0), 1e-30))
                            scvfFluxes[scvfIndexInInside[localScvfIdx]] = 0;
                        // cubes
                        else if (dim == 1 || geomType.isCube())
                        {
                            const auto fIdx = scvfIndexInInside[localScvfIdx];
                            const auto fIdxOpposite = fIdx%2 ? fIdx-1 : fIdx+1;
                            scvfFluxes[fIdx] = -scvfFluxes[fIdxOpposite];
                        }
                        // simplices
                        else if (geomType.isSimplex())
                            scvfFluxes[scvfIndexInInside[localScvfIdx]] = 0;
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
                DUNE_THROW(Dune::NotImplemented, "velocity output for cell-centered and prism/pyramid");

            Velocity scvVelocity(0);
            jacobianT2.mtv(refVelocity, scvVelocity);

            scvVelocity /= geometry.integrationElement(localPos);

            int eIdxGlobal = fvGridGeometry_.elementMapper().index(element);

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
    // The following SFINAE enable_if usage allows compilation, even if only a
    //
    // boundaryTypes(const Element&, const scv&)
    //
    // is provided in the problem file. In that case, the compiler cannot detect
    // (without additional measures like "using...") the signature
    //
    // boundaryTypes(const Element&, const scvf&)
    //
    // in the problem base class. Therefore, calls to this method trigger a
    // compiler error. However, that call is needed for calculating velocities
    // if the cell-centered discretization is used. By proceeding as in the
    // following lines, that call will only be compiled if cell-centered
    // actually is used.
    template <bool enable = isBox, typename std::enable_if_t<!enable, int> = 0>
    BoundaryTypes problemBoundaryTypes_(const Element& element, const SubControlVolumeFace& scvf) const
    { return problem_.boundaryTypes(element, scvf); }

    //! we should never call this method for box models
    template <bool enable = isBox, typename std::enable_if_t<enable, int> = 0>
    BoundaryTypes problemBoundaryTypes_(const Element& element, const SubControlVolumeFace& scvf) const
    { return BoundaryTypes(); }

    const Problem& problem_;
    const FVGridGeometry& fvGridGeometry_;
    const GridVariables& gridVariables_;

    bool enableOutput_;
    std::vector<int> cellNum_;
};

} // end namespace Dumux

#endif
