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
 * \brief Velocity output for implicit (porous media) models
 */
#ifndef DUMUX_IMPLICIT_VELOCITYOUTPUT_HH
#define DUMUX_IMPLICIT_VELOCITYOUTPUT_HH

#include <dumux/implicit/properties.hh>
#include "problem.hh"
#include <unordered_map>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

namespace Dumux
{

//At the moment this property is defined in the individual models -> should be changed
namespace Properties
{
    NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
}

/*!
 * \brief Velocity output for implicit (porous media) models
 */
template<class TypeTag>
class ImplicitVelocityOutput
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    ImplicitVelocityOutput(const Problem& problem)
    : problem_(problem)
    {
        // check, if velocity output can be used (works only for cubes so far)
        velocityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity);
        if (velocityOutput_)
        {
            if (isBox)
            {
                cellNum_.assign(problem_.gridView().size(dofCodim), 0);

                for (const auto& element : elements(problem_.gridView()))
                {
                    FVElementGeometry fvGeometry;
                    fvGeometry.update(problem_.gridView(), element);

                    for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                    {
                        int vIdxGlobal = problem_.vertexMapper().subIndex(element, scvIdx, dofCodim);
                        cellNum_[vIdxGlobal] += 1;
                    }
                }
            }
        }
    }

    bool enableOutput()
    {
        return velocityOutput_;
    }

    // The following SFINAE enable_if usage allows compilation, even if only a
    //
    // boundaryTypes(BoundaryTypes&, const Vertex&)
    //
    // is provided in the problem file. In that case, the compiler cannot detect
    // (without additional measures like "using...") the signature
    //
    // boundaryTypes(BoundaryTypes&, const Intersection&)
    //
    // in the problem base class. Therefore, calls to this method trigger a
    // compiler error. However, that call is needed for calculating velocities
    // if the cell-centered discretization is used. By proceeding as in the
    // following lines, that call will only be compiled if cell-centered
    // actually is used.
    template <class T = TypeTag>
    void problemBoundaryTypes(BoundaryTypes& bcTypes,
                              const typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox), Intersection>::type& intersection) const
    {
        problem_.boundaryTypes(bcTypes, intersection);
    }
    template <class T = TypeTag>
    void problemBoundaryTypes(BoundaryTypes& bcTypes,
                              const typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox), Intersection>::type& intersection) const
    {}

    template<class VelocityVector>
    void calculateVelocity(VelocityVector& velocity,
                           const ElementVolumeVariables& elemVolVars,
                           const FVElementGeometry& fvGeometry,
                           const Element& element,
                           int phaseIdx)
    {
        if (velocityOutput_)
        {
            const auto geometry = element.geometry();

            Dune::GeometryType geomType = geometry.type();
            const ReferenceElement &referenceElement
                = ReferenceElements::general(geomType);

            const Dune::FieldVector<Scalar, dim>& localPos
                = referenceElement.position(0, 0);

            // get the transposed Jacobian of the element mapping
            const typename Element::Geometry::JacobianTransposed jacobianT2 =
                geometry.jacobianTransposed(localPos);

            if (isBox)
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > ScvVelocities;
                ScvVelocities scvVelocities(fvGeometry.numScv);
                scvVelocities = 0;

                for (int fIdx = 0; fIdx < fvGeometry.numScvf; fIdx++)
                {
                    // local position of integration point
                    const Dune::FieldVector<Scalar, dim>& localPosIP = fvGeometry.subContVolFace[fIdx].ipLocal;

                    // Transformation of the global normal vector to normal vector in the reference element
                    const typename Element::Geometry::JacobianTransposed jacobianT1 =
                        geometry.jacobianTransposed(localPosIP);

                    FluxVariables fluxVars(problem_,
                                           element,
                                           fvGeometry,
                                           fIdx,
                                           elemVolVars);

                    const GlobalPosition globalNormal = fluxVars.face().normal;

                    GlobalPosition localNormal(0);
                    jacobianT1.mv(globalNormal, localNormal);
                    localNormal /= localNormal.two_norm();

                    // area of the subcontrolvolume face in the reference element
                    Scalar localArea = scvfReferenceArea_(geomType, fIdx);

                    // get the volume flux divided by the area of the
                    // subcontrolvolume face in the reference element
                    Scalar flux = fluxVars.volumeFlux(phaseIdx) / localArea;

                    // transform the volume flux into a velocity vector
                    GlobalPosition tmpVelocity = localNormal;
                    tmpVelocity *= flux;

                    scvVelocities[fluxVars.face().i] += tmpVelocity;
                    scvVelocities[fluxVars.face().j] += tmpVelocity;
                }

                // transform vertex velocities from local to global coordinates
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int vIdxGlobal = problem_.vertexMapper().subIndex(element, scvIdx, dofCodim);

                    // calculate the subcontrolvolume velocity by the Piola transformation
                    Dune::FieldVector<CoordScalar, dimWorld> scvVelocity(0);

                    jacobianT2.mtv(scvVelocities[scvIdx], scvVelocity);
                    scvVelocity /= geometry.integrationElement(localPos)*cellNum_[vIdxGlobal];
                    // add up the wetting phase subcontrolvolume velocities for each vertex
                    velocity[vIdxGlobal] += scvVelocity;
                }
            }
            else
            {
                std::vector<Scalar> scvfFluxes(element.subEntities(1), 0);

                int fIdxInner = 0;
                for (const auto& intersection : intersections(problem_.gridView(), element))
                {
                    int fIdx = intersection.indexInInside();

                    if (intersection.neighbor())
                    {
                        FluxVariables fluxVars(problem_,
                                               element,
                                               fvGeometry,
                                               fIdxInner,
                                               elemVolVars);

                        scvfFluxes[fIdx] += fluxVars.volumeFlux(phaseIdx);

                        fIdxInner++;
                    }
                    else if (intersection.boundary())
                    {
                        FluxVariables fluxVars(problem_,
                                               element,
                                               fvGeometry,
                                               fIdx,
                                               elemVolVars,true);

                        scvfFluxes[fIdx] = fluxVars.volumeFlux(phaseIdx);
                    }
                }

                // Correct boundary fluxes in case of Neumann conditions.
                // They are simply set to an average of the inner fluxes. In
                // this general setting, it would be very difficult to
                // calculate correct phase, i.e., volume, fluxes from arbitrary
                // Neumann conditions.
                if (element.hasBoundaryIntersections())
                {
                    for (const auto& intersection : intersections(problem_.gridView(), element))
                    {
                        if (intersection.boundary())
                        {
                            BoundaryTypes bcTypes;
                            problemBoundaryTypes(bcTypes, intersection);

                            if (bcTypes.hasNeumann())
                            {
                                int fIdx = intersection.indexInInside();

                                // cubes
                                if (dim == 1 || geomType.isCube()){
                                    int fIdxOpposite = fIdx%2 ? fIdx-1 : fIdx+1;
                                    scvfFluxes[fIdx] = -scvfFluxes[fIdxOpposite];
                                }
                                // simplices
                                else if (geomType.isSimplex()) {
                                    scvfFluxes[fIdx] = 0;
                                }
                            }
                        }
                    }
                }

                Dune::FieldVector < Scalar, dim > refVelocity;
                // cubes
                if (dim == 1 || geomType.isCube()){
                    for (int i = 0; i < dim; i++)
                    {
                        refVelocity[i] = 0.5 * (scvfFluxes[2*i + 1] - scvfFluxes[2*i]);
                    }
                }
                // simplices: Raviart-Thomas-0 interpolation evaluated at the cell center
                else if (geomType.isSimplex()) {
                    for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                    {
                        refVelocity[dimIdx] = -scvfFluxes[dim - 1 - dimIdx];
                         for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                         {
                             refVelocity[dimIdx] += scvfFluxes[fIdx]/(dim + 1);
                         }
                    }
                }
                // 3D prism and pyramids
                else {
                    DUNE_THROW(Dune::NotImplemented,
                               "velocity output for cell-centered and prism/pyramid");
                }

                Dune::FieldVector<Scalar, dimWorld> scvVelocity(0);
                jacobianT2.mtv(refVelocity, scvVelocity);

                scvVelocity /= geometry.integrationElement(localPos);

                int eIdxGlobal = problem_.elementMapper().index(element);

                velocity[eIdxGlobal]= scvVelocity;
            }
        } // velocity output
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

protected:
    const Problem& problem_;
    bool velocityOutput_;
    std::vector<int> cellNum_;
};

}
#endif
