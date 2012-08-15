/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
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

#ifndef DUMUX_MPFAL2PFABOUNDVELOCITY2P_HH
#define DUMUX_MPFAL2PFABOUNDVELOCITY2P_HH

#include "fvmpfal2pfaboundpressure2p.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Markus Wolff
 */

namespace Dumux
{

template<class TypeTag> class FVMPFAL2PFABoundVelocity2P: public FVMPFAL2PFABoundPressure2P<TypeTag>
{
    typedef FVMPFAL2PFABoundPressure2P<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElementContainer;
    typedef Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    typedef typename GET_PROP_TYPE(TypeTag, GridTypeIndices) GridTypeIndices;

    typedef typename ParentType::InteractionVolume InteractionVolume;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    enum
    {
        globalCorner = 2,
        globalEdge = 3,
        neumannNeumann = 0,
        dirichletDirichlet = 1,
        dirichletNeumann = 2,
        neumannDirichlet = 3
    };

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    FVMPFAL2PFABoundVelocity2P(Problem& problem) :
            ParentType(problem), problem_(problem), gravity_(problem.gravity())
    {
        ElementIterator element = problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(*element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(*element));
        fluidState.setTemperature(problem_.temperature(*element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);
    }

    void calculateVelocity();

    void initialize(bool solveTwice = true)
    {
        ParentType::initialize(solveTwice);

        calculateVelocity();

        return;
    }

    void update()
    {
        ParentType::update();

        calculateVelocity();

        return;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);

        Dune::BlockVector < DimVector > &velocityWetting = *(writer.template allocateManagedBuffer<
                Scalar, dim>(problem_.gridView().size(0)));
        Dune::BlockVector < DimVector > &velocityNonwetting =
                *(writer.template allocateManagedBuffer<Scalar, dim>(problem_.gridView().size(0)));

        // compute update vector
        ElementIterator eItEnd = problem_.gridView().template end<0>();
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // cell index
            int globalIdx = problem_.variables().index(*eIt);

            CellData& cellData = problem_.variables().cellData(globalIdx);

            Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
            Dune::FieldVector < Scalar, 2 * dim > fluxNW(0);
            // run through all intersections with neighbors and boundary
            IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                int isIndex = isIt->indexInInside();

                fluxW[isIndex] = isIt->geometry().volume()
                        * (isIt->centerUnitOuterNormal() * cellData.fluxData().velocity(wPhaseIdx, isIndex));
                fluxNW[isIndex] =
                        isIt->geometry().volume()
                                * (isIt->centerUnitOuterNormal()
                                        * cellData.fluxData().velocity(nPhaseIdx, isIndex));
            }

            DimVector refVelocity(0);
            refVelocity[0] = 0.5 * (fluxW[1] - fluxW[0]);
            refVelocity[1] = 0.5 * (fluxW[3] - fluxW[2]);

            const DimVector& localPos =
                    ReferenceElementContainer::general(eIt->geometry().type()).position(0, 0);

            // get the transposed Jacobian of the element mapping
            const DimMatrix& jacobianT = eIt->geometry().jacobianTransposed(localPos);

            // calculate the element velocity by the Piola transformation
            DimVector elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= eIt->geometry().integrationElement(localPos);

            velocityWetting[globalIdx] = elementVelocity;

            refVelocity = 0;
            refVelocity[0] = 0.5 * (fluxNW[1] - fluxNW[0]);
            refVelocity[1] = 0.5 * (fluxNW[3] - fluxNW[2]);

            // calculate the element velocity by the Piola transformation
            elementVelocity = 0;
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= eIt->geometry().integrationElement(localPos);

            velocityNonwetting[globalIdx] = elementVelocity;
        }

        writer.attachCellData(velocityWetting, "wetting-velocity", dim);
        writer.attachCellData(velocityNonwetting, "non-wetting-velocity", dim);

        return;
    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    static const Scalar threshold_ = 1e-15;
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation); //!< gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
};
// end of template

template<class TypeTag>
void FVMPFAL2PFABoundVelocity2P<TypeTag>::calculateVelocity()
{
    // run through all elements
    VertexIterator vItEnd = problem_.gridView().template end<dim>();
    for (VertexIterator vIt = problem_.gridView().template begin<dim>(); vIt != vItEnd; ++vIt)
    {
        int globalVertIdx = problem_.variables().index(*vIt);

        InteractionVolume& interactionVolume = this->interactionVolumes_[globalVertIdx];

        if (interactionVolume.isInnerVolume())
        {
            ElementPointer & elementPointer1 = interactionVolume.getSubVolumeElement(0);
            ElementPointer & elementPointer2 = interactionVolume.getSubVolumeElement(1);
            ElementPointer & elementPointer3 = interactionVolume.getSubVolumeElement(2);
            ElementPointer & elementPointer4 = interactionVolume.getSubVolumeElement(3);

            // cell index
            int globalIdx1 = problem_.variables().index(*elementPointer1);
            int globalIdx2 = problem_.variables().index(*elementPointer2);
            int globalIdx3 = problem_.variables().index(*elementPointer3);
            int globalIdx4 = problem_.variables().index(*elementPointer4);

            //get the cell Data
            CellData& cellData1 = problem_.variables().cellData(globalIdx1);
            CellData& cellData2 = problem_.variables().cellData(globalIdx2);
            CellData& cellData3 = problem_.variables().cellData(globalIdx3);
            CellData& cellData4 = problem_.variables().cellData(globalIdx4);

            // get global coordinate of cell centers
            const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
            const GlobalPosition& globalPos2 = elementPointer2->geometry().center();
            const GlobalPosition& globalPos3 = elementPointer3->geometry().center();
            const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

            // get pressure values
            Dune::FieldVector < Scalar, 2 * dim > pW(0);
            Dune::FieldVector < Scalar, 2 * dim > pN(0);

            pW[0] = cellData1.pressure(wPhaseIdx);
            pW[1] = cellData2.pressure(wPhaseIdx);
            pW[2] = cellData3.pressure(wPhaseIdx);
            pW[3] = cellData4.pressure(wPhaseIdx);

            pN[0] = cellData1.pressure(nPhaseIdx);
            pN[1] = cellData2.pressure(nPhaseIdx);
            pN[2] = cellData3.pressure(nPhaseIdx);
            pN[3] = cellData4.pressure(nPhaseIdx);

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda1(cellData1.mobility(wPhaseIdx));
            lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda2(cellData2.mobility(wPhaseIdx));
            lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda3(cellData3.mobility(wPhaseIdx));
            lambda3[nPhaseIdx] = cellData3.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal3 = lambda3[wPhaseIdx] + lambda3[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda4(cellData4.mobility(wPhaseIdx));
            lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

            std::vector < DimVector > lambda(2 * dim);

            lambda[0][0] = lambdaTotal1;
            lambda[0][1] = lambdaTotal1;
            lambda[1][0] = lambdaTotal2;
            lambda[1][1] = lambdaTotal2;
            lambda[2][0] = lambdaTotal3;
            lambda[2][1] = lambdaTotal3;
            lambda[3][0] = lambdaTotal4;
            lambda[3][1] = lambdaTotal4;

            Scalar potentialW12 = 0;
            Scalar potentialW14 = 0;
            Scalar potentialW32 = 0;
            Scalar potentialW34 = 0;

            Scalar potentialNW12 = 0;
            Scalar potentialNW14 = 0;
            Scalar potentialNW32 = 0;
            Scalar potentialNW34 = 0;

            //flux vector
            Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
            Dune::FieldVector < Scalar, 2 * dim > fluxNW(0);

            Dune::FieldMatrix < Scalar, dim, 2 * dim - dim + 1 > T(0);
            DimVector Tu(0);
            Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);

            bool rightTriangle = this->calculateTransmissibility(T, interactionVolume, lambda, 0, 1, 2, 3);

            if (rightTriangle)
            {
                u[0] = pW[1];
                u[1] = pW[2];
                u[2] = pW[0];

                T.mv(u, Tu);

                fluxW[0] = Tu[1];
                potentialW12 = Tu[1];

                u[0] = pN[1];
                u[1] = pN[2];
                u[2] = pN[0];

                T.mv(u, Tu);

                fluxNW[0] = Tu[1];
                potentialNW12 = Tu[1];
            }
            else
            {
                u[0] = pW[0];
                u[1] = pW[3];
                u[2] = pW[1];

                T.mv(u, Tu);

                fluxW[0] = Tu[1];
                potentialW12 = Tu[1];

                u[0] = pN[0];
                u[1] = pN[3];
                u[2] = pN[1];

                T.mv(u, Tu);

                fluxNW[0] = Tu[1];
                potentialNW12 = Tu[1];
            }

            rightTriangle = this->calculateTransmissibility(T, interactionVolume, lambda, 1, 2, 3, 0);

            if (rightTriangle)
            {
                u[0] = pW[2];
                u[1] = pW[3];
                u[2] = pW[1];

                T.mv(u, Tu);

                fluxW[1] = Tu[1];
                potentialW32 = -Tu[1];

                u[0] = pN[2];
                u[1] = pN[3];
                u[2] = pN[1];

                T.mv(u, Tu);

                fluxNW[1] = Tu[1];
                potentialNW32 = -Tu[1];
            }
            else
            {
                u[0] = pW[1];
                u[1] = pW[0];
                u[2] = pW[2];

                T.mv(u, Tu);

                fluxW[1] = Tu[1];
                potentialW32 = -Tu[1];

                u[0] = pN[1];
                u[1] = pN[0];
                u[2] = pN[2];

                T.mv(u, Tu);

                fluxNW[1] = Tu[1];
                potentialNW32 = -Tu[1];
            }

            rightTriangle = this->calculateTransmissibility(T, interactionVolume, lambda, 2, 3, 0, 1);

            if (rightTriangle)
            {
                u[0] = pW[3];
                u[1] = pW[0];
                u[2] = pW[2];

                T.mv(u, Tu);

                fluxW[2] = Tu[1];
                potentialW34 = Tu[1];

                u[0] = pN[3];
                u[1] = pN[0];
                u[2] = pN[2];

                T.mv(u, Tu);

                fluxNW[2] = Tu[1];
                potentialNW34 = Tu[1];
            }
            else
            {
                u[0] = pW[2];
                u[1] = pW[1];
                u[2] = pW[3];

                T.mv(u, Tu);

                fluxW[2] = Tu[1];
                potentialW34 = Tu[1];

                u[0] = pN[2];
                u[1] = pN[1];
                u[2] = pN[3];

                T.mv(u, Tu);

                fluxNW[2] = Tu[1];
                potentialNW34 = Tu[1];
            }

            rightTriangle = this->calculateTransmissibility(T, interactionVolume, lambda, 3, 0, 1, 2);

            if (rightTriangle)
            {
                u[0] = pW[0];
                u[1] = pW[1];
                u[2] = pW[3];

                T.mv(u, Tu);

                fluxW[3] = Tu[1];
                potentialW14 = -Tu[1];

                u[0] = pN[0];
                u[1] = pN[1];
                u[2] = pN[3];

                T.mv(u, Tu);

                fluxNW[3] = Tu[1];
                potentialNW14 = -Tu[1];
            }
            else
            {
                u[0] = pW[3];
                u[1] = pW[2];
                u[2] = pW[0];

                T.mv(u, Tu);

                fluxW[3] = Tu[1];
                potentialW14 = -Tu[1];

                u[0] = pN[3];
                u[1] = pN[2];
                u[2] = pN[0];

                T.mv(u, Tu);

                fluxNW[3] = Tu[1];
                potentialNW14 = -Tu[1];
            }

            //store potentials for further calculations (saturation, ...)
            cellData1.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialW12);
            cellData1.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialNW12);
            cellData1.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 1), potentialW14);
            cellData1.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 1), potentialNW14);
            cellData2.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 0), -potentialW32);
            cellData2.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 0), -potentialNW32);
            cellData2.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialW12);
            cellData2.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialNW12);
            cellData3.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 0), potentialW34);
            cellData3.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 0), potentialNW34);
            cellData3.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 1), potentialW32);
            cellData3.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 1), potentialNW32);
            cellData4.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -potentialW14);
            cellData4.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -potentialNW14);
            cellData4.fluxData().addPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialW34);
            cellData4.fluxData().addPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialNW34);

            //compute mobilities of face 1
            Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
            lambda12Upw[wPhaseIdx] = (potentialW12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
            lambda12Upw[nPhaseIdx] = (potentialNW12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

            //compute mobilities of face 4
            Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
            lambda14Upw[wPhaseIdx] = (potentialW14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
            lambda14Upw[nPhaseIdx] = (potentialNW14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

            //compute mobilities of face 2
            Dune::FieldVector<Scalar, numPhases> lambda32Upw(0.0);
            lambda32Upw[wPhaseIdx] = (potentialW32 >= 0) ? lambda3[wPhaseIdx] : lambda2[wPhaseIdx];
            lambda32Upw[nPhaseIdx] = (potentialNW32 >= 0) ? lambda3[nPhaseIdx] : lambda2[nPhaseIdx];

            //compute mobilities of face 3
            Dune::FieldVector<Scalar, numPhases> lambda34Upw(0.0);
            lambda34Upw[wPhaseIdx] = (potentialW34 >= 0) ? lambda3[wPhaseIdx] : lambda4[wPhaseIdx];
            lambda34Upw[nPhaseIdx] = (potentialNW34 >= 0) ? lambda3[nPhaseIdx] : lambda4[nPhaseIdx];

            for (int i = 0; i < numPhases; i++)
            {
                // evaluate parts of velocity --> always take the normal for which the flux is calculated!
                DimVector vel12 = interactionVolume.getNormal(0, 0);
                DimVector vel14 = interactionVolume.getNormal(3, 0);
                DimVector vel23 = interactionVolume.getNormal(1, 0);
                DimVector vel21 = interactionVolume.getNormal(0, 0);
                DimVector vel34 = interactionVolume.getNormal(2, 0);
                DimVector vel32 = interactionVolume.getNormal(1, 0);
                DimVector vel41 = interactionVolume.getNormal(3, 0);
                DimVector vel43 = interactionVolume.getNormal(2, 0);

                Dune::FieldVector<Scalar, 2 * dim> flux(0);
                switch (i)
                {
                case wPhaseIdx:
                {
                    flux = fluxW;
                    break;
                }
                case nPhaseIdx:
                {
                    flux = fluxNW;
                    break;
                }
                }

                vel12 *= flux[0] / (2 * interactionVolume.getFaceArea(0, 0)); //divide by 2 because the flux is related to the half face!
                vel14 *= flux[3] / (2 * interactionVolume.getFaceArea(0, 1));
                vel23 *= flux[1] / (2 * interactionVolume.getFaceArea(1, 0));
                vel21 *= flux[0] / (2 * interactionVolume.getFaceArea(1, 1));
                vel34 *= flux[2] / (2 * interactionVolume.getFaceArea(2, 0));
                vel32 *= flux[1] / (2 * interactionVolume.getFaceArea(2, 1));
                vel41 *= flux[3] / (2 * interactionVolume.getFaceArea(3, 0));
                vel43 *= flux[2] / (2 * interactionVolume.getFaceArea(3, 1));

                Scalar lambdaT12 = lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx];
                Scalar lambdaT14 = lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx];
                Scalar lambdaT32 = lambda32Upw[wPhaseIdx] + lambda32Upw[nPhaseIdx];
                Scalar lambdaT34 = lambda34Upw[wPhaseIdx] + lambda34Upw[nPhaseIdx];
                Scalar fracFlow12 = (lambdaT12 > threshold_) ? lambda12Upw[i] / (lambdaT12) : 0.0;
                Scalar fracFlow14 = (lambdaT14 > threshold_) ? lambda14Upw[i] / (lambdaT14) : 0.0;
                Scalar fracFlow32 = (lambdaT32 > threshold_) ? lambda32Upw[i] / (lambdaT32) : 0.0;
                Scalar fracFlow34 = (lambdaT34 > threshold_) ? lambda34Upw[i] / (lambdaT34) : 0.0;

                vel12 *= fracFlow12;
                vel14 *= fracFlow14;
                vel23 *= fracFlow32;
                vel21 *= fracFlow12;
                vel34 *= fracFlow34;
                vel32 *= fracFlow32;
                vel41 *= fracFlow14;
                vel43 *= fracFlow34;

                if (this->innerBoundaryVolumeFaces_[globalIdx1][interactionVolume.getIndexOnElement(0, 0)])
                {
                    vel12 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx1][interactionVolume.getIndexOnElement(0, 1)])
                {
                    vel14 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx2][interactionVolume.getIndexOnElement(1, 0)])
                {
                    vel23 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx2][interactionVolume.getIndexOnElement(1, 1)])
                {
                    vel21 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx3][interactionVolume.getIndexOnElement(2, 0)])
                {
                    vel34 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx3][interactionVolume.getIndexOnElement(2, 1)])
                {
                    vel32 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx4][interactionVolume.getIndexOnElement(3, 0)])
                {
                    vel41 *= 2;
                }
                if (this->innerBoundaryVolumeFaces_[globalIdx4][interactionVolume.getIndexOnElement(3, 1)])
                {
                    vel43 *= 2;
                }

                //store velocities
                cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 0), vel12);
                cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 1), vel14);
                cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 0), vel23);
                cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 1), vel21);
                cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 0), vel34);
                cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 1), vel32);
                cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 0), vel41);
                cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 1), vel43);
            }
            //set velocity marker
            cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 0));
            cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 1));
            cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 0));
            cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 1));
            cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 0));
            cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 1));
            cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 0));
            cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 1));
        }

        // at least one face on boundary!
        else
        {
            for (int elemIdx = 0; elemIdx < 2 * dim; elemIdx++)
            {
                bool isOutside = false;
                for (int faceIdx = 0; faceIdx < dim; faceIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, faceIdx);
                    if (interactionVolume.isOutsideFace(intVolFaceIdx))
                    {
                        isOutside = true;
                        break;
                    }
                }
                if (isOutside)
                {
                    continue;
                }

                ElementPointer & elementPointer = interactionVolume.getSubVolumeElement(elemIdx);

                // get global coordinate of cell centers
                const GlobalPosition& globalPos = elementPointer->geometry().center();

                // cell index
                int globalIdx = problem_.variables().index(*elementPointer);

                //get the cell Data
                CellData& cellData = problem_.variables().cellData(globalIdx);

                //permeability vector at boundary
                DimMatrix permeability(problem_.spatialParams().intrinsicPermeability(*elementPointer));

                //get mobilities of the phases
                Dune::FieldVector<Scalar, numPhases> lambda(cellData.mobility(wPhaseIdx));
                lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

                Scalar pressW = cellData.pressure(wPhaseIdx);
                Scalar pressNW = cellData.pressure(nPhaseIdx);

                for (int faceIdx = 0; faceIdx < dim; faceIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, faceIdx);

                    if (interactionVolume.isBoundaryFace(intVolFaceIdx))
                    {
                        if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(pressEqIdx))
                        {
                            int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, faceIdx);

                            const ReferenceElement& referenceElement = ReferenceElementContainer::general(
                                    elementPointer->geometry().type());

                            const LocalPosition& localPos = referenceElement.position(boundaryFaceIdx, 1);

                            const GlobalPosition& globalPosFace = elementPointer->geometry().global(localPos);

                            DimVector distVec(globalPosFace - globalPos);
                            Scalar dist = distVec.two_norm();
                            DimVector unitDistVec(distVec);
                            unitDistVec /= dist;

                            // get pc and lambda at the boundary
                            Scalar satWBound = cellData.saturation(wPhaseIdx);
                            //check boundary sat at face 1
                            if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(satEqIdx))
                            {
                                Scalar satBound = interactionVolume.getDirichletValues(intVolFaceIdx)[saturationIdx];
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satWBound = satBound;
                                      break;
                                }
                                case Sn:
                                {
                                    satWBound = 1 - satBound;
                                    break;
                                }
                                }

                            }

                            Scalar pcBound = MaterialLaw::pC(
                                    problem_.spatialParams().materialLawParams(*elementPointer), satWBound);

                            Scalar gravityDiffBound = (problem_.bboxMax() - globalPosFace) * problem_.gravity()
                                    * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                            pcBound += gravityDiffBound;

                            Dune::FieldVector<Scalar, numPhases> lambdaBound(
                                    MaterialLaw::krw(problem_.spatialParams().materialLawParams(*elementPointer),
                                            satWBound));
                            lambdaBound[nPhaseIdx] = MaterialLaw::krn(
                                    problem_.spatialParams().materialLawParams(*elementPointer), satWBound);
                            lambdaBound[wPhaseIdx] /= viscosity_[wPhaseIdx];
                            lambdaBound[nPhaseIdx] /= viscosity_[nPhaseIdx];

                            Scalar gdeltaZ = (problem_.bboxMax()-globalPos) * gravity_;
                            Scalar potentialBoundW = interactionVolume.getDirichletValues(intVolFaceIdx)[pressureIdx] + density_[wPhaseIdx]*gdeltaZ;
                            Scalar potentialBoundNW = potentialBoundW;

                            //calculate potential gradients
                            switch (pressureType_)
                            {
                            case pw:
                            {
                                potentialBoundNW += pcBound;
                                break;
                            }
                            case pn:
                            {
                                //calculate potential gradients
                                potentialBoundW -= pcBound;
                                break;
                            }
                            }


                            Scalar potentialW = (pressW - potentialBoundW) / dist;
                            Scalar  potentialNW = (pressNW - potentialBoundNW) / dist;

                            //store potentials for further calculations (saturation, ...)
                            cellData.fluxData().addPotential(wPhaseIdx, boundaryFaceIdx, potentialW);
                            cellData.fluxData().addPotential(nPhaseIdx, boundaryFaceIdx, potentialNW);

                            //calculated phase velocities from advective velocities -> capillary pressure velocity already added in pressure part!
                            DimVector velocityW(0);
                            DimVector velocityNW(0);

                            // calculate capillary pressure gradient
                            DimVector pressGradient = unitDistVec;
                            pressGradient *= (pressW - potentialBoundW) / dist;
                            permeability.mv(pressGradient, velocityW);

                            pressGradient = unitDistVec;
                            pressGradient *= (pressNW - potentialBoundNW) / dist;
                            permeability.mv(pressGradient, velocityNW);

                            velocityW *= (potentialW >= 0.) ? lambda[wPhaseIdx] : lambdaBound[wPhaseIdx];
                            velocityNW *= (potentialNW >= 0.) ? lambda[nPhaseIdx] : lambdaBound[nPhaseIdx];

                            //velocity is calculated from two vertices of one intersection!
                            velocityW *= 0.5;
                            velocityNW *= 0.5;

                            //store velocities
                                velocityW += cellData.fluxData().velocity(wPhaseIdx, boundaryFaceIdx);
                                velocityNW += cellData.fluxData().velocity(nPhaseIdx, boundaryFaceIdx);
                                cellData.fluxData().setVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                                cellData.fluxData().setVelocity(nPhaseIdx, boundaryFaceIdx, velocityNW);
                                cellData.fluxData().setVelocityMarker(boundaryFaceIdx);

                        }
                        else if (interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx))
                        {
                            int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, faceIdx);

                            const ReferenceElement& referenceElement = ReferenceElementContainer::general(
                                    elementPointer->geometry().type());

                            const LocalPosition& localPos = referenceElement.position(boundaryFaceIdx, 1);

                            const GlobalPosition& globalPosFace = elementPointer->geometry().global(localPos);

                            DimVector distVec(globalPosFace - globalPos);
                            Scalar dist = distVec.two_norm();
                            DimVector unitDistVec(distVec);
                            unitDistVec /= dist;

                            // get neumann boundary value
                            PrimaryVariables boundValues(interactionVolume.getNeumannValues(intVolFaceIdx));

                            boundValues[wPhaseIdx] /= density_[wPhaseIdx];
                            boundValues[nPhaseIdx] /= density_[nPhaseIdx];

                            DimVector velocityW(unitDistVec);
                            DimVector velocityNW(unitDistVec);

                            velocityW *= boundValues[wPhaseIdx] / (2 * interactionVolume.getFaceArea(elemIdx, faceIdx));
                            velocityNW *= boundValues[nPhaseIdx]
                                    / (2 * interactionVolume.getFaceArea(elemIdx, faceIdx));

                            //store potentials for further calculations (saturation, ...)
                            cellData.fluxData().addPotential(wPhaseIdx, boundaryFaceIdx, boundValues[wPhaseIdx]);
                            cellData.fluxData().addPotential(nPhaseIdx, boundaryFaceIdx, boundValues[nPhaseIdx]);

                            //store velocities
                            velocityW += cellData.fluxData().velocity(wPhaseIdx, boundaryFaceIdx);
                            velocityNW += cellData.fluxData().velocity(nPhaseIdx, boundaryFaceIdx);
                            cellData.fluxData().setVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                            cellData.fluxData().setVelocity(nPhaseIdx, boundaryFaceIdx, velocityNW);
                            cellData.fluxData().setVelocityMarker(boundaryFaceIdx);
                        }
                        else
                        {
                            DUNE_THROW(Dune::NotImplemented,
                                    "No valid boundary condition type defined for pressure equation!");
                        }
                    }
                }
            }

        } // end boundaries

    } // end vertex iterator

//    }
    return;
} // end method calcTotalVelocity

}
// end of Dune namespace
#endif
