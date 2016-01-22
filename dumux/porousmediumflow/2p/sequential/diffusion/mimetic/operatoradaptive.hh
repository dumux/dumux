// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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

#ifndef DUMUX_MIMETICOPERATOR2PADAPTIVE_HH
#define DUMUX_MIMETICOPERATOR2PADAPTIVE_HH

/*!
 * \file
 *
 * \brief An assembler for the Jacobian matrix based on mimetic FD.
 */

#include "croperatoradaptive.hh"
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/mimetic/properties.hh>

namespace Dumux
{
/*!
 * \ingroup Mimetic2P
 * @brief Levelwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - Scalar The field type used in the elements of the stiffness matrix
*/
template<class TypeTag>
class MimeticOperatorAssemblerTwoPAdaptive : public CROperatorAssemblerTwoPAdaptive<TypeTag>
{
    typedef CROperatorAssemblerTwoPAdaptive<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum
    {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    typedef typename GET_PROP_TYPE(TypeTag, LocalStiffness) LocalStiffness;

    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation),
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        saturationIdx = Indices::saturationIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef Dune::FieldVector<Scalar, dimWorld> FieldVector;
public:

    MimeticOperatorAssemblerTwoPAdaptive (const GridView& gridView)
    : ParentType(gridView)
    {}

    template<class Vector>
    void calculatePressure(LocalStiffness& loc, Vector& u, Problem& problem)
    {
        Dune::DynamicVector<Scalar> velocityW(2*dim);
        Dune::DynamicVector<Scalar> velocityNw(2*dim);
        Dune::DynamicVector<Scalar> pressTraceW(2*dim);
        Dune::DynamicVector<Scalar> pressTraceNw(2*dim);

        const auto firstElement = *problem.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem.referencePressure(firstElement));
        fluidState.setPressure(nPhaseIdx, problem.referencePressure(firstElement));
        fluidState.setTemperature(problem.temperature(firstElement));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        Scalar densityDiff = FluidSystem::density(fluidState, nPhaseIdx) - FluidSystem::density(fluidState, wPhaseIdx);
        Scalar viscosityW = FluidSystem::viscosity(fluidState, wPhaseIdx);
        Scalar viscosityNw = FluidSystem::viscosity(fluidState, nPhaseIdx);

        //reset velocity
        for (int i = 0; i < problem.gridView().size(0); i++)
        {
            problem.variables().cellData(i).fluxData().resetVelocity();
        }

        // run over all level elements
        for (const auto& element : Dune::elements(this->gridView_))
        {
            int eIdxGlobal = problem.variables().index(element);

            unsigned int numFaces = this->intersectionMapper_.size(eIdxGlobal);

            // get local to global id map and pressure traces
            velocityW.resize(numFaces);
            velocityNw.resize(numFaces);
            pressTraceW.resize(numFaces);
            pressTraceNw.resize(numFaces);

            CellData& cellData = problem.variables().cellData(eIdxGlobal);
            FieldVector globalPos = element.geometry().center();

            int intersectionIdx = -1;
            // get local to global id map and pressure traces
            for (const auto& intersection : Dune::intersections(problem.gridView(), element))
            {
                ++intersectionIdx;

                int fIdxGlobal = this->intersectionMapper_.subIndex(element, intersectionIdx);

                Scalar pcPotFace = (problem.bBoxMax() - intersection.geometry().center()) * problem.gravity() * densityDiff;

                switch (pressureType)
                {
                case pw:
                {
                    pressTraceW[intersectionIdx] = u[fIdxGlobal];
                    pressTraceNw[intersectionIdx] = u[fIdxGlobal] + pcPotFace;
                    break;
                }
                case pn:
                {
                    pressTraceNw[intersectionIdx] = u[fIdxGlobal];
                    pressTraceW[intersectionIdx] = u[fIdxGlobal] - pcPotFace;

                    break;
                }
                }
            }

            switch (pressureType)
            {
            case pw:
            {
                Scalar potW = loc.constructPressure(element, pressTraceW);
                Scalar gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * densityDiff;
                Scalar potNw = potW + gravPot;

                cellData.setPotential(wPhaseIdx, potW);
                cellData.setPotential(nPhaseIdx, potNw);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, wPhaseIdx);

                cellData.setPressure(wPhaseIdx, potW - gravPot);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, nPhaseIdx);

                cellData.setPressure(nPhaseIdx, potNw - gravPot);

                break;
            }
            case pn:
            {
                Scalar potNw = loc.constructPressure(element, pressTraceNw);
                Scalar  gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * densityDiff;
                Scalar potW = potNw - gravPot;

                cellData.setPotential(nPhaseIdx, potNw);
                cellData.setPotential(wPhaseIdx, potW);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, wPhaseIdx);

                cellData.setPressure(wPhaseIdx, potW - gravPot);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, nPhaseIdx);

                cellData.setPressure(nPhaseIdx, potNw - gravPot);

                break;
            }
            }

            //velocity reconstruction: !!! The velocity which is not reconstructed from the primary
            //pressure variable can be slightly wrong and not conservative!!!!
            // -> Should not be used for transport!!
            loc.constructVelocity(element, velocityW, pressTraceW, cellData.potential(wPhaseIdx));
            loc.constructVelocity(element, velocityNw, pressTraceNw, cellData.potential(nPhaseIdx));

            intersectionIdx = -1;
            for (const auto& intersection : Dune::intersections(problem.gridView(), element))
            {
                ++intersectionIdx;
                int idxInInside = intersection.indexInInside();

                cellData.fluxData().addUpwindPotential(wPhaseIdx, idxInInside, velocityW[intersectionIdx]);
                cellData.fluxData().addUpwindPotential(nPhaseIdx, idxInInside, velocityNw[intersectionIdx]);

                Scalar mobilityW = 0;
                Scalar mobilityNw = 0;

                if (intersection.neighbor())
                {
                    int neighborIdx = problem.variables().index(intersection.outside());

                    CellData& cellDataNeighbor = problem.variables().cellData(neighborIdx);

                    mobilityW =
                            (velocityW[intersectionIdx] >= 0.) ? cellData.mobility(wPhaseIdx) :
                                    cellDataNeighbor.mobility(wPhaseIdx);
                    mobilityNw =
                            (velocityNw[intersectionIdx] >= 0.) ? cellData.mobility(nPhaseIdx) :
                                    cellDataNeighbor.mobility(nPhaseIdx);

                    if (velocityW[intersectionIdx] >= 0.)
                    {
                        FieldVector velocity(intersection.centerUnitOuterNormal());
                        velocity *= mobilityW/(mobilityW+mobilityNw) * velocityW[intersectionIdx];
                        cellData.fluxData().addVelocity(wPhaseIdx, idxInInside, velocity);
                        cellDataNeighbor.fluxData().addVelocity(wPhaseIdx, intersection.indexInOutside(), velocity);
                    }
                    if (velocityNw[intersectionIdx] >= 0.)
                    {
                        FieldVector velocity(intersection.centerUnitOuterNormal());
                        velocity *= mobilityNw/(mobilityW+mobilityNw) * velocityNw[intersectionIdx];
                        cellData.fluxData().addVelocity(nPhaseIdx, idxInInside, velocity);
                        cellDataNeighbor.fluxData().addVelocity(nPhaseIdx, intersection.indexInOutside(), velocity);
                    }

                    cellData.fluxData().setVelocityMarker(idxInInside);
                }
                else
                {
                    BoundaryTypes bctype;
                    problem.boundaryTypes(bctype, intersection);
                    if (bctype.isDirichlet(satEqIdx))
                    {
                        PrimaryVariables boundValues(0.0);
                        problem.dirichlet(boundValues, intersection);

                        if (velocityW[intersectionIdx] >= 0.)
                        {
                            mobilityW = cellData.mobility(wPhaseIdx);
                        }
                        else
                        {
                            mobilityW = MaterialLaw::krw(problem.spatialParams().materialLawParams(element),
                                boundValues[saturationIdx]) / viscosityW;
                        }

                        if (velocityNw[intersectionIdx] >= 0.)
                        {
                            mobilityNw = cellData.mobility(nPhaseIdx);
                        }
                        else
                        {
                            mobilityNw = MaterialLaw::krn(problem.spatialParams().materialLawParams(element),
                                boundValues[saturationIdx]) / viscosityNw;
                        }
                    }
                    else
                    {
                        mobilityW = cellData.mobility(wPhaseIdx);
                        mobilityNw = cellData.mobility(nPhaseIdx);
                    }

                    FieldVector velocity(intersection.centerUnitOuterNormal());
                    velocity *= mobilityW/(mobilityW+mobilityNw) * velocityW[intersectionIdx];
                    cellData.fluxData().addVelocity(wPhaseIdx, idxInInside, velocity);


                    velocity = intersection.centerUnitOuterNormal();
                    velocity *= mobilityNw/(mobilityW+mobilityNw) * velocityNw[intersectionIdx];
                    cellData.fluxData().addVelocity(nPhaseIdx, idxInInside, velocity);
                    cellData.fluxData().setVelocityMarker(idxInInside);
                }
            }
        }
    }
};
}
#endif
