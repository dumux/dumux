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
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 */
#ifndef DUMUX_IMPLICIT_VELOCITYOUTPUT_HH
#define DUMUX_IMPLICIT_VELOCITYOUTPUT_HH

#include "implicitproperties.hh"
#include <unordered_map>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/istl/bvector.hh>

namespace Dumux
{

//At the moment this property is defined in the individual models -> should be changed
namespace Properties
{
    NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
}

template<class TypeTag>
class ImplicitVelocityOutput
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;

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
        if (velocityOutput_ && isBox)
        {
            cellNum_.assign(problem_.gridView().size(dofCodim), 0);
        }

        ElementIterator eIt = problem_.gridView().template begin<0>();
        ElementIterator eEndIt = problem_.gridView().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            if (eIt->geometry().type().isCube() == false)
            {
                velocityOutput_ = false;
            }

            if (velocityOutput_ && isBox)
            {
                FVElementGeometry fvGeometry;
                fvGeometry.update(problem_.gridView(), *eIt);

                // transform vertex velocities from local to global coordinates
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                    int vIdxGlobal = problem_.vertexMapper().subIndex(*eIt, scvIdx, dofCodim);
#else
                    int vIdxGlobal = problem_.vertexMapper().map(*eIt, scvIdx, dofCodim);
#endif

                    cellNum_[vIdxGlobal] += 1;
                }
            }
        }

        if (velocityOutput_ != GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity))
            std::cout << "ATTENTION: Velocity output only works for cubes and is set to false for simplices\n";
    }

    bool enableOutput()
    {
        return velocityOutput_;
    }

    template<class VelocityVector>
    void calculateVelocity(VelocityVector& velocity,
                           const ElementVolumeVariables& elemVolVars,
                           const FVElementGeometry& fvGeometry,
                           const Element& element, 
                           int phaseIdx)
    {
        if (velocityOutput_)
        {
            // calculate vertex velocities
            GlobalPosition tmpVelocity(0.0);

            const Dune::FieldVector<Scalar, dim>& localPos =
                ReferenceElements::general(element.geometry().type()).position(0, 0);

            // get the transposed Jacobian of the element mapping
            const typename Element::Geometry::JacobianTransposed jacobianT2 =
                element.geometry().jacobianTransposed(localPos);

            if (isBox)
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > SCVVelocities;
                SCVVelocities scvVelocities(pow(2,dim));
                scvVelocities = 0;

                for (int fIdx = 0; fIdx < fvGeometry.numScvf; fIdx++)
                {
                    // local position of integration point
                    const Dune::FieldVector<Scalar, dim>& localPosIP = fvGeometry.subContVolFace[fIdx].ipLocal;

                    // Transformation of the global normal vector to normal vector in the reference element
                    const typename Element::Geometry::JacobianTransposed jacobianT1 =
                        element.geometry().jacobianTransposed(localPosIP);

                    FluxVariables fluxVars(problem_,
                                           element,
                                           fvGeometry,
                                           fIdx,
                                           elemVolVars);

                    const GlobalPosition globalNormal = fluxVars.face().normal;

                    GlobalPosition localNormal(0);
                    jacobianT1.mv(globalNormal, localNormal);
                    // note only works for cubes
                    const Scalar localArea = pow(2,-(dim-1));

                    localNormal /= localNormal.two_norm();

                    // Get the Darcy velocities. The Darcy velocities are divided by the area of the subcontrolvolume
                    // face in the reference element.
                    Scalar flux = fluxVars.volumeFlux(phaseIdx) / localArea;

                    // transform the normal Darcy velocity into a vector
                    tmpVelocity = localNormal;
                    tmpVelocity *= flux;

                    scvVelocities[fluxVars.face().i] += tmpVelocity;
                    scvVelocities[fluxVars.face().j] += tmpVelocity;
                }

                // transform vertex velocities from local to global coordinates
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                    int vIdxGlobal = problem_.vertexMapper().subIndex(element, scvIdx, dofCodim);
#else
                    int vIdxGlobal = problem_.vertexMapper().map(element, scvIdx, dofCodim);
#endif
                    // calculate the subcontrolvolume velocity by the Piola transformation
                    Dune::FieldVector<CoordScalar, dim> scvVelocity(0);

                    jacobianT2.mtv(scvVelocities[scvIdx], scvVelocity);
                    scvVelocity /= element.geometry().integrationElement(localPos)*cellNum_[vIdxGlobal];
                    // add up the wetting phase subcontrolvolume velocities for each vertex
                    velocity[vIdxGlobal] += scvVelocity;
                }
            }
            else
            {
                Dune::FieldVector<Scalar, 2*dim> scvVelocities(0.0);

                int fIdxInner = 0;
                IntersectionIterator isEndIt = problem_.gridView().iend(element);
                for (IntersectionIterator isIt = problem_.gridView().ibegin(element); 
                     isIt != isEndIt; ++isIt)
                {
                    int fIdx = isIt->indexInInside();

                    if (isIt->neighbor())
                    {
                        FluxVariables fluxVars(problem_,
                                               element,
                                               fvGeometry,
                                               fIdxInner,
                                               elemVolVars);

                        Scalar flux = fluxVars.volumeFlux(phaseIdx);
                        scvVelocities[fIdx] += flux;

                        fIdxInner++;
                    }
                    else if (isIt->boundary())
                    {
                        FluxVariables fluxVars(problem_,
                                               element,
                                               fvGeometry,
                                               fIdx,
                                               elemVolVars,true);

                        Scalar flux = fluxVars.volumeFlux(phaseIdx);
                        scvVelocities[fIdx] = flux;
                    }
                }

                Dune::FieldVector < Scalar, dim > refVelocity(0);
                for (int i = 0; i < dim; i++)
                    refVelocity[i] = 0.5 * (scvVelocities[2*i + 1] - scvVelocities[2*i]);

                Dune::FieldVector<Scalar, dimWorld> scvVelocity(0);
                jacobianT2.mtv(refVelocity, scvVelocity);

                scvVelocity /= element.geometry().integrationElement(localPos);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                int eIdxGlobal = problem_.elementMapper().index(element);
#else
                int eIdxGlobal = problem_.elementMapper().map(element);
#endif

                velocity[eIdxGlobal]= scvVelocity;
            }
        } // velocity output
    }

protected:
    const Problem& problem_;
    bool velocityOutput_;
    std::vector<int> cellNum_;
};

}
#endif
