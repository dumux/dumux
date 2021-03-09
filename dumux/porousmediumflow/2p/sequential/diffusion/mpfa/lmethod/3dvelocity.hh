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
 * \ingroup SequentialTwoPModel
 * \brief 2-d velocity calculation using a 3-d MPFA L-method.
 */
#ifndef DUMUX_FVMPFAL2PVELOCITY2P_HH
#define DUMUX_FVMPFAL2PVELOCITY2P_HH

#include <dune/grid/common/gridenums.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>
#include "3dtransmissibilitycalculator.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class for calculating 3-d velocities from cell-wise constant pressure values.
 *
 * Calculates phase velocities or total velocity from a known pressure field applying a finite volume discretization and a MPFA L-method.
 * At Dirichlet boundaries a two-point flux approximation is used.
 * The pressure has to be given as piecewise constant cell values.
 * The velocities are calculated as
 *
 *\f[ \boldsymbol v_\alpha = - \lambda_\alpha \boldsymbol K \textbf{grad}\, \Phi_\alpha, \f]
 * and,
 * \f[ \boldsymbol v_t = \boldsymbol v_w + \boldsymbol v_n,\f]
 *
 * where \f$ \Phi_\alpha \f$ denotes the potential of phase \f$ \alpha \f$, \f$ \boldsymbol K \f$ the intrinsic permeability,
 * and \f$ \lambda_\alpha \f$ a phase mobility.
 *
 * Remark1: only for 3-D Hexahedrons of quadrilaterals!
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag> class FvMpfaL3dVelocity2p
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Grid = typename GridView::Grid;
    using IndexSet = typename GridView::IndexSet;

    using Geometry = typename Element::Geometry;
    using JacobianTransposed = typename Geometry::JacobianTransposed ;

    using GridTypeIndices = GetPropType<TypeTag, Properties::GridTypeIndices>;

    using InteractionVolume = GetPropType<TypeTag, Properties::MPFAInteractionVolume>;
    using InteractionVolumeContainer = GetPropType<TypeTag, Properties::MPFAInteractionVolumeContainer>;
    using TransmissibilityCalculator = FvMpfaL3dTransmissibilityCalculator<TypeTag>;
    using TransmissibilityType = typename TransmissibilityCalculator::TransmissibilityType;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        sw = Indices::saturationW,
        sn = Indices::saturationNw,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        vt = Indices::velocityTotal
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressureEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
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
    enum
    {
        innerEdgeFace = 2, innerSideFace = 1
    };

    using GlobalPosition = typename Geometry::GlobalCoordinate;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    /*!
     * \brief Constructs a FvMpfaL3dVelocity2p object
     *
     * \param problem A problem class object
     */
    FvMpfaL3dVelocity2p(Problem& problem) :
        problem_(problem), gravity_(problem.gravity())
    {
        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        vtkOutputLevel_ = getParam<int>("Vtk.OutputLevel");
    }

    //! Calculate velocities for flux faces of an interaction volume
    void calculateInnerInteractionVolumeVelocity(InteractionVolume& interactionVolume,
            CellData & cellData1,  CellData & cellData2, CellData & cellData3, CellData & cellData4,
            CellData & cellData5, CellData & cellData6, CellData & cellData7, CellData & cellData8,
            InteractionVolumeContainer& interactionVolumes,
            TransmissibilityCalculator& transmissibilityCalculator, int fIdx = -1);
    void calculateBoundaryInteractionVolumeVelocity(InteractionVolume& interactionVolume,
                                                    CellData& cellData, int elemIdx);

    //! Initializes the velocity model
    void initialize(bool solveTwice = true)
    {
        const auto element = *problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
        fluidState.setTemperature(problem_.temperature(element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);

        return;
    }

    /*!
     * \brief Adds velocity output to the output file
     *
     * Adds the phase velocities or a total velocity (depending on the formulation) to the output.
     * If the VtkOutputLevel is equal to zero (default) only primary variables are written,
     * if it is larger than zero also secondary variables are written.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        if (vtkOutputLevel_ > 0)
        {
            Dune::BlockVector < DimVector > &velocityWetting = *(writer.template allocateManagedBuffer<Scalar,
                                                                 dim>(problem_.gridView().size(0)));
            Dune::BlockVector < DimVector > &velocityNonwetting = *(writer.template allocateManagedBuffer<Scalar,
                                                                    dim>(problem_.gridView().size(0)));

            // compute update vector
            for (const auto& element : elements(problem_.gridView()))
            {
                // cell index
                int eIdxGlobal = problem_.variables().index(element);

                CellData & cellData = problem_.variables().cellData(eIdxGlobal);

                Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
                Dune::FieldVector < Scalar, 2 * dim > fluxNw(0);

                // run through all intersections with neighbors and boundary
                for (const auto& intersection : intersections(problem_.gridView(), element))
                {
                    int isIndex = intersection.indexInInside();

                    fluxW[isIndex] += intersection.geometry().volume()
                                * (intersection.centerUnitOuterNormal() * cellData.fluxData().velocity(wPhaseIdx, isIndex));
                    fluxNw[isIndex] += intersection.geometry().volume()
                                * (intersection.centerUnitOuterNormal() * cellData.fluxData().velocity(nPhaseIdx, isIndex));
                }

                DimVector refVelocity(0);
                refVelocity[0] = 0.5 * (fluxW[1] - fluxW[0]);
                refVelocity[1] = 0.5 * (fluxW[3] - fluxW[2]);
                refVelocity[2] = 0.5 * (fluxW[5] - fluxW[4]);

                const DimVector& localPos = referenceElement(element).position(0, 0);

                // get the transposed Jacobian of the element mapping
                const JacobianTransposed jacobianT = element.geometry().jacobianTransposed(localPos);

                // calculate the element velocity by the Piola transformation
                DimVector elementVelocity(0);
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= element.geometry().integrationElement(localPos);

                velocityWetting[eIdxGlobal] = elementVelocity;

                refVelocity = 0;
                refVelocity[0] = 0.5 * (fluxNw[1] - fluxNw[0]);
                refVelocity[1] = 0.5 * (fluxNw[3] - fluxNw[2]);
                refVelocity[2] = 0.5 * (fluxNw[5] - fluxNw[4]);

                // calculate the element velocity by the Piola transformation
                elementVelocity = 0;
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= element.geometry().integrationElement(localPos);

                velocityNonwetting[eIdxGlobal] = elementVelocity;
            }

            writer.attachCellData(velocityWetting, "wetting-velocity", dim);
            writer.attachCellData(velocityNonwetting, "nonwetting-velocity", dim);
        }

        return;
    }

private:
    Problem& problem_;
    const Dune::FieldVector<Scalar, dimWorld>& gravity_; //!< vector including the gravity constant

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;

    static constexpr Scalar threshold_ = 1e-15;
    //! gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int velocityType_ = getPropValue<TypeTag, Properties::VelocityFormulation>();
    //! gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();
    //! gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
};
// end of template

/*!
 * \brief Calculates the velocities at the flux faces of an interation volume around a vertex which is not a boundary vertex.
 *
 *  Calculates the velocities at the flux faces of an interation volume around a vertex
 *  which is not a boundary vertex and adds them to the face velocity vectors in the <tt>CellData</tt> objects.
 *
 * \param interactionVolume An <tt>InteractionVolume</tt> object including the information for calculating the MPFA transmissibilities
 * \param cellData1  <tt>CellData</tt> object of an IMPES model for sub-volume 1
 * \param cellData2  <tt>CellData</tt> object of an IMPES model for sub-volume 2
 * \param cellData3  <tt>CellData</tt> object of an IMPES model for sub-volume 3
 * \param cellData4  <tt>CellData</tt> object of an IMPES model for sub-volume 4
 * \param cellData5  <tt>CellData</tt> object of an IMPES model for sub-volume 5
 * \param cellData6  <tt>CellData</tt> object of an IMPES model for sub-volume 6
 * \param cellData7  <tt>CellData</tt> object of an IMPES model for sub-volume 7
 * \param cellData8  <tt>CellData</tt> object of an IMPES model for sub-volume 8
 * \param interactionVolumes Container including the interaction volume information for the complete grid
 * \param transmissibilityCalculator Object including the methods for calculating the transmissibilities
 * \param fIdx Index of the flux face for which the velocity has to be calculated. If no face index is given, <tt>fIdx</tt> = -1
 * and velocities for all flux faces in the interaction volume are calculated!
 */
template<class TypeTag>
void FvMpfaL3dVelocity2p<TypeTag>::calculateInnerInteractionVolumeVelocity(InteractionVolume& interactionVolume,
        CellData & cellData1,  CellData & cellData2, CellData & cellData3, CellData & cellData4,
        CellData & cellData5, CellData & cellData6, CellData & cellData7, CellData & cellData8,
        InteractionVolumeContainer& interactionVolumes, TransmissibilityCalculator& transmissibilityCalculator, int fIdx)
        {
    auto element1 = interactionVolume.getSubVolumeElement(0);
    auto element2 = interactionVolume.getSubVolumeElement(1);
    auto element3 = interactionVolume.getSubVolumeElement(2);
    auto element4 = interactionVolume.getSubVolumeElement(3);
    auto element5 = interactionVolume.getSubVolumeElement(4);
    auto element6 = interactionVolume.getSubVolumeElement(5);
    auto element7 = interactionVolume.getSubVolumeElement(6);
    auto element8 = interactionVolume.getSubVolumeElement(7);

    // cell index
    int eIdxGlobal1 = problem_.variables().index(element1);
    int eIdxGlobal2 = problem_.variables().index(element2);
    int eIdxGlobal3 = problem_.variables().index(element3);
    int eIdxGlobal4 = problem_.variables().index(element4);
    int eIdxGlobal5 = problem_.variables().index(element5);
    int eIdxGlobal6 = problem_.variables().index(element6);
    int eIdxGlobal7 = problem_.variables().index(element7);
    int eIdxGlobal8 = problem_.variables().index(element8);

    // pressures flux calculation
    Dune::FieldVector<Scalar, 8> potW(0);
    Dune::FieldVector<Scalar, 8> potNw(0);

    potW[0] = cellData1.potential(wPhaseIdx);
    potW[1] = cellData2.potential(wPhaseIdx);
    potW[2] = cellData3.potential(wPhaseIdx);
    potW[3] = cellData4.potential(wPhaseIdx);
    potW[4] = cellData5.potential(wPhaseIdx);
    potW[5] = cellData6.potential(wPhaseIdx);
    potW[6] = cellData7.potential(wPhaseIdx);
    potW[7] = cellData8.potential(wPhaseIdx);

    potNw[0] = cellData1.potential(nPhaseIdx);
    potNw[1] = cellData2.potential(nPhaseIdx);
    potNw[2] = cellData3.potential(nPhaseIdx);
    potNw[3] = cellData4.potential(nPhaseIdx);
    potNw[4] = cellData5.potential(nPhaseIdx);
    potNw[5] = cellData6.potential(nPhaseIdx);
    potNw[6] = cellData7.potential(nPhaseIdx);
    potNw[7] = cellData8.potential(nPhaseIdx);

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda1(cellData1.mobility(wPhaseIdx));
    lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

    // compute total mobility of cell 1
    Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda2(cellData2.mobility(wPhaseIdx));
    lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

    // compute total mobility of cell 2
    Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda3(cellData3.mobility(wPhaseIdx));
    lambda3[nPhaseIdx] = cellData3.mobility(nPhaseIdx);

    // compute total mobility of cell 3
    Scalar lambdaTotal3 = lambda3[wPhaseIdx] + lambda3[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda4(cellData4.mobility(wPhaseIdx));
    lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

    // compute total mobility of cell 4
    Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda5(cellData5.mobility(wPhaseIdx));
    lambda5[nPhaseIdx] = cellData5.mobility(nPhaseIdx);

    // compute total mobility of cell 5
    Scalar lambdaTotal5 = lambda5[wPhaseIdx] + lambda5[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda6(cellData6.mobility(wPhaseIdx));
    lambda6[nPhaseIdx] = cellData6.mobility(nPhaseIdx);

    // compute total mobility of cell 6
    Scalar lambdaTotal6 = lambda6[wPhaseIdx] + lambda6[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda7(cellData7.mobility(wPhaseIdx));
    lambda7[nPhaseIdx] = cellData7.mobility(nPhaseIdx);

    // compute total mobility of cell 7
    Scalar lambdaTotal7 = lambda7[wPhaseIdx] + lambda7[nPhaseIdx];

    // get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda8(cellData8.mobility(wPhaseIdx));
    lambda8[nPhaseIdx] = cellData8.mobility(nPhaseIdx);

    // compute total mobility of cell 8
    Scalar lambdaTotal8 = lambda8[wPhaseIdx] + lambda8[nPhaseIdx];

    std::vector<DimVector> lambda(8);
    lambda[0][0] = lambdaTotal1;
    lambda[0][1] = lambdaTotal1;
    lambda[0][2] = lambdaTotal1;
    lambda[1][0] = lambdaTotal2;
    lambda[1][1] = lambdaTotal2;
    lambda[1][2] = lambdaTotal2;
    lambda[2][0] = lambdaTotal3;
    lambda[2][1] = lambdaTotal3;
    lambda[2][2] = lambdaTotal3;
    lambda[3][0] = lambdaTotal4;
    lambda[3][1] = lambdaTotal4;
    lambda[3][2] = lambdaTotal4;
    lambda[4][0] = lambdaTotal5;
    lambda[4][1] = lambdaTotal5;
    lambda[4][2] = lambdaTotal5;
    lambda[5][0] = lambdaTotal6;
    lambda[5][1] = lambdaTotal6;
    lambda[5][2] = lambdaTotal6;
    lambda[6][0] = lambdaTotal7;
    lambda[6][1] = lambdaTotal7;
    lambda[6][2] = lambdaTotal7;
    lambda[7][0] = lambdaTotal8;
    lambda[7][1] = lambdaTotal8;
    lambda[7][2] = lambdaTotal8;

    Scalar potentialDiffW0 = 0;
    Scalar potentialDiffW1 = 0;
    Scalar potentialDiffW2 = 0;
    Scalar potentialDiffW3 = 0;
    Scalar potentialDiffW4 = 0;
    Scalar potentialDiffW5 = 0;
    Scalar potentialDiffW6 = 0;
    Scalar potentialDiffW7 = 0;
    Scalar potentialDiffW8 = 0;
    Scalar potentialDiffW9 = 0;
    Scalar potentialDiffW10 = 0;
    Scalar potentialDiffW11 = 0;

    Scalar potentialDiffNw0 = 0;
    Scalar potentialDiffNw1 = 0;
    Scalar potentialDiffNw2 = 0;
    Scalar potentialDiffNw3 = 0;
    Scalar potentialDiffNw4 = 0;
    Scalar potentialDiffNw5 = 0;
    Scalar potentialDiffNw6 = 0;
    Scalar potentialDiffNw7 = 0;
    Scalar potentialDiffNw8 = 0;
    Scalar potentialDiffNw9 = 0;
    Scalar potentialDiffNw10 = 0;
    Scalar potentialDiffNw11 = 0;

    //flux vector
    Dune::FieldVector<Scalar, 12> fluxW(0);
    Dune::FieldVector<Scalar, 12> fluxNw(0);

    DimVector Tu(0);
    Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);
    TransmissibilityType T(0);

    if (fIdx < 0 || fIdx == 0)
    {
        // calculate the flux through the subvolumeface 1 (subVolumeFaceIdx = 0)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 0, 1, 2, 3,
                4, 5);

        if (caseL == 1)
        {
            u[0] = potW[0];
            u[1] = potW[1];
            u[2] = potW[2];
            u[3] = potW[4];

            T.mv(u, Tu);

            fluxW[0] = Tu[0];
            potentialDiffW0 = Tu[0];

            u[0] = potNw[0];
            u[1] = potNw[1];
            u[2] = potNw[2];
            u[3] = potNw[4];

            T.mv(u, Tu);

            fluxNw[0] = Tu[0];
            potentialDiffNw0 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[0];
            u[1] = potW[1];
            u[2] = potW[3];
            u[3] = potW[5];

            T.mv(u, Tu);

            fluxW[0] = Tu[0];
            potentialDiffW0 = Tu[0];

            u[0] = potNw[0];
            u[1] = potNw[1];
            u[2] = potNw[3];
            u[3] = potNw[5];

            T.mv(u, Tu);

            fluxNw[0] = Tu[0];
            potentialDiffNw0 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[0];
            u[1] = potW[1];
            u[2] = potW[3];
            u[3] = potW[4];

            T.mv(u, Tu);

            fluxW[0] = Tu[0];
            potentialDiffW0 = Tu[0];

            u[0] = potNw[0];
            u[1] = potNw[1];
            u[2] = potNw[3];
            u[3] = potNw[4];

            T.mv(u, Tu);

            fluxNw[0] = Tu[0];
            potentialDiffNw0 = Tu[0];
        }
        else
        {
            u[0] = potW[0];
            u[1] = potW[1];
            u[2] = potW[2];
            u[3] = potW[5];

            T.mv(u, Tu);

            fluxW[0] = Tu[0];
            potentialDiffW0 = Tu[0];

            u[0] = potNw[0];
            u[1] = potNw[1];
            u[2] = potNw[2];
            u[3] = potNw[5];

            T.mv(u, Tu);

            fluxNw[0] = Tu[0];
            potentialDiffNw0 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 1)
    {
        // calculate the flux through the subvolumeface 2 (subVolumeFaceIdx = 1)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 1, 3, 0, 2, 5,
                7);

        if (caseL == 1)
        {
            u[0] = potW[1];
            u[1] = potW[3];
            u[2] = potW[0];
            u[3] = potW[5];

            T.mv(u, Tu);

            fluxW[1] = Tu[0];
            potentialDiffW1 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[3];
            u[2] = potNw[0];
            u[3] = potNw[5];

            T.mv(u, Tu);

            fluxNw[1] = Tu[0];
            potentialDiffNw1 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[1];
            u[1] = potW[3];
            u[2] = potW[2];
            u[3] = potW[7];

            T.mv(u, Tu);

            fluxW[1] = Tu[0];
            potentialDiffW1 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[3];
            u[2] = potNw[2];
            u[3] = potNw[7];

            T.mv(u, Tu);

            fluxNw[1] = Tu[0];
            potentialDiffNw1 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[1];
            u[1] = potW[3];
            u[2] = potW[2];
            u[3] = potW[5];

            T.mv(u, Tu);

            fluxW[1] = Tu[0];
            potentialDiffW1 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[3];
            u[2] = potNw[2];
            u[3] = potNw[5];

            T.mv(u, Tu);

            fluxNw[1] = Tu[0];
            potentialDiffNw1 = Tu[0];
        }
        else
        {
            u[0] = potW[1];
            u[1] = potW[3];
            u[2] = potW[0];
            u[3] = potW[7];

            T.mv(u, Tu);

            fluxW[1] = Tu[0];
            potentialDiffW1 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[3];
            u[2] = potNw[0];
            u[3] = potNw[7];

            T.mv(u, Tu);

            fluxNw[1] = Tu[0];
            potentialDiffNw1 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 2)
    {
        // calculate the flux through the subvolumeface 3 (subVolumeFaceIdx = 2)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 3, 2, 1, 0, 7,
                6);

        if (caseL == 1)
        {
            u[0] = potW[3];
            u[1] = potW[2];
            u[2] = potW[1];
            u[3] = potW[7];

            T.mv(u, Tu);

            fluxW[2] = Tu[0];
            potentialDiffW2 = Tu[0];

            u[0] = potNw[3];
            u[1] = potNw[2];
            u[2] = potNw[1];
            u[3] = potNw[7];

            T.mv(u, Tu);

            fluxNw[2] = Tu[0];
            potentialDiffNw2 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[3];
            u[1] = potW[2];
            u[2] = potW[0];
            u[3] = potW[6];

            T.mv(u, Tu);

            fluxW[2] = Tu[0];
            potentialDiffW2 = Tu[0];

            u[0] = potNw[3];
            u[1] = potNw[2];
            u[2] = potNw[0];
            u[3] = potNw[6];

            T.mv(u, Tu);

            fluxNw[2] = Tu[0];
            potentialDiffNw2 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[3];
            u[1] = potW[2];
            u[2] = potW[0];
            u[3] = potW[7];

            T.mv(u, Tu);

            fluxW[2] = Tu[0];
            potentialDiffW2 = Tu[0];

            u[0] = potNw[3];
            u[1] = potNw[2];
            u[2] = potNw[0];
            u[3] = potNw[7];

            T.mv(u, Tu);

            fluxNw[2] = Tu[0];
            potentialDiffNw2 = Tu[0];
        }
        else
        {
            u[0] = potW[3];
            u[1] = potW[2];
            u[2] = potW[1];
            u[3] = potW[6];

            T.mv(u, Tu);

            fluxW[2] = Tu[0];
            potentialDiffW2 = Tu[0];

            u[0] = potNw[3];
            u[1] = potNw[2];
            u[2] = potNw[1];
            u[3] = potNw[6];

            T.mv(u, Tu);

            fluxNw[2] = Tu[0];
            potentialDiffNw2 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 3)
    {
        // calculate the flux through the subvolumeface 4 (subVolumeFaceIdx = 3)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 2, 0, 3, 1, 6,
                4);

        if (caseL == 1)
        {
            u[0] = potW[2];
            u[1] = potW[0];
            u[2] = potW[3];
            u[3] = potW[6];

            T.mv(u, Tu);

            fluxW[3] = Tu[0];
            potentialDiffW3 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[0];
            u[2] = potNw[3];
            u[3] = potNw[6];

            T.mv(u, Tu);

            fluxNw[3] = Tu[0];
            potentialDiffNw3 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[2];
            u[1] = potW[0];
            u[2] = potW[1];
            u[3] = potW[4];

            T.mv(u, Tu);

            fluxW[3] = Tu[0];
            potentialDiffW3 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[0];
            u[2] = potNw[1];
            u[3] = potNw[4];

            T.mv(u, Tu);

            fluxNw[3] = Tu[0];
            potentialDiffNw3 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[2];
            u[1] = potW[0];
            u[2] = potW[1];
            u[3] = potW[6];

            T.mv(u, Tu);

            fluxW[3] = Tu[0];
            potentialDiffW3 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[0];
            u[2] = potNw[1];
            u[3] = potNw[6];

            T.mv(u, Tu);

            fluxNw[3] = Tu[0];
            potentialDiffNw3 = Tu[0];
        }
        else
        {
            u[0] = potW[2];
            u[1] = potW[0];
            u[2] = potW[3];
            u[3] = potW[4];

            T.mv(u, Tu);

            fluxW[3] = Tu[0];
            potentialDiffW3 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[0];
            u[2] = potNw[3];
            u[3] = potNw[4];

            T.mv(u, Tu);

            fluxNw[3] = Tu[0];
            potentialDiffNw3 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 4)
    {
        // calculate the flux through the subvolumeface 5 (subVolumeFaceIdx = 4)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 5, 4, 7, 6, 1,
                0);

        if (caseL == 1)
        {
            u[0] = potW[5];
            u[1] = potW[4];
            u[2] = potW[7];
            u[3] = potW[1];

            T.mv(u, Tu);

            fluxW[4] = Tu[0];
            potentialDiffW4 = Tu[0];

            u[0] = potNw[5];
            u[1] = potNw[4];
            u[2] = potNw[7];
            u[3] = potNw[1];

            T.mv(u, Tu);

            fluxNw[4] = Tu[0];
            potentialDiffNw4 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[5];
            u[1] = potW[4];
            u[2] = potW[6];
            u[3] = potW[0];

            T.mv(u, Tu);

            fluxW[4] = Tu[0];
            potentialDiffW4 = Tu[0];

            u[0] = potNw[5];
            u[1] = potNw[4];
            u[2] = potNw[6];
            u[3] = potNw[0];

            T.mv(u, Tu);

            fluxNw[4] = Tu[0];
            potentialDiffNw4 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[5];
            u[1] = potW[4];
            u[2] = potW[6];
            u[3] = potW[1];

            T.mv(u, Tu);

            fluxW[4] = Tu[0];
            potentialDiffW4 = Tu[0];

            u[0] = potNw[5];
            u[1] = potNw[4];
            u[2] = potNw[6];
            u[3] = potNw[1];

            T.mv(u, Tu);

            fluxNw[4] = Tu[0];
            potentialDiffNw4 = Tu[0];
        }
        else
        {
            u[0] = potW[5];
            u[1] = potW[4];
            u[2] = potW[7];
            u[3] = potW[0];

            T.mv(u, Tu);

            fluxW[4] = Tu[0];
            potentialDiffW4 = Tu[0];

            u[0] = potNw[5];
            u[1] = potNw[4];
            u[2] = potNw[7];
            u[3] = potNw[0];

            T.mv(u, Tu);

            fluxNw[4] = Tu[0];
            potentialDiffNw4 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 5)
    {
        // calculate the flux through the subvolumeface 6 (subVolumeFaceIdx = 5)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 7, 5, 6, 4, 3,
                1);

        if (caseL == 1)
        {
            u[0] = potW[7];
            u[1] = potW[5];
            u[2] = potW[6];
            u[3] = potW[3];

            T.mv(u, Tu);

            fluxW[5] = Tu[0];
            potentialDiffW5 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[5];
            u[2] = potNw[6];
            u[3] = potNw[3];

            T.mv(u, Tu);

            fluxNw[5] = Tu[0];
            potentialDiffNw5 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[7];
            u[1] = potW[5];
            u[2] = potW[4];
            u[3] = potW[1];

            T.mv(u, Tu);

            fluxW[5] = Tu[0];
            potentialDiffW5 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[5];
            u[2] = potNw[4];
            u[3] = potNw[1];

            T.mv(u, Tu);

            fluxNw[5] = Tu[0];
            potentialDiffNw5 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[7];
            u[1] = potW[5];
            u[2] = potW[4];
            u[3] = potW[3];

            T.mv(u, Tu);

            fluxW[5] = Tu[0];
            potentialDiffW5 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[5];
            u[2] = potNw[4];
            u[3] = potNw[3];

            T.mv(u, Tu);

            fluxNw[5] = Tu[0];
            potentialDiffNw5 = Tu[0];
        }
        else
        {
            u[0] = potW[7];
            u[1] = potW[5];
            u[2] = potW[6];
            u[3] = potW[1];

            T.mv(u, Tu);

            fluxW[5] = Tu[0];
            potentialDiffW5 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[5];
            u[2] = potNw[6];
            u[3] = potNw[1];

            T.mv(u, Tu);

            fluxNw[5] = Tu[0];
            potentialDiffNw5 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 6)
    {
        // calculate the flux through the subvolumeface 7 (subVolumeFaceIdx = 6)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 6, 7, 4, 5, 2,
                3);

        if (caseL == 1)
        {
            u[0] = potW[6];
            u[1] = potW[7];
            u[2] = potW[4];
            u[3] = potW[2];

            T.mv(u, Tu);

            fluxW[6] = Tu[0];
            potentialDiffW6 = Tu[0];

            u[0] = potNw[6];
            u[1] = potNw[7];
            u[2] = potNw[4];
            u[3] = potNw[2];

            T.mv(u, Tu);

            fluxNw[6] = Tu[0];
            potentialDiffNw6 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[6];
            u[1] = potW[7];
            u[2] = potW[5];
            u[3] = potW[3];

            T.mv(u, Tu);

            fluxW[6] = Tu[0];
            potentialDiffW6 = Tu[0];

            u[0] = potNw[6];
            u[1] = potNw[7];
            u[2] = potNw[5];
            u[3] = potNw[3];

            T.mv(u, Tu);

            fluxNw[6] = Tu[0];
            potentialDiffNw6 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[6];
            u[1] = potW[7];
            u[2] = potW[5];
            u[3] = potW[2];

            T.mv(u, Tu);

            fluxW[6] = Tu[0];
            potentialDiffW6 = Tu[0];

            u[0] = potNw[6];
            u[1] = potNw[7];
            u[2] = potNw[5];
            u[3] = potNw[2];

            T.mv(u, Tu);

            fluxNw[6] = Tu[0];
            potentialDiffNw6 = Tu[0];
        }
        else
        {
            u[0] = potW[6];
            u[1] = potW[7];
            u[2] = potW[4];
            u[3] = potW[3];

            T.mv(u, Tu);

            fluxW[6] = Tu[0];
            potentialDiffW6 = Tu[0];

            u[0] = potNw[6];
            u[1] = potNw[7];
            u[2] = potNw[4];
            u[3] = potNw[3];

            T.mv(u, Tu);

            fluxNw[6] = Tu[0];
            potentialDiffNw6 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 7)
    {
        // calculate the flux through the subvolumeface 8 (subVolumeFaceIdx = 7)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 4, 6, 5, 7, 0,
                2);

        if (caseL == 1)
        {
            u[0] = potW[4];
            u[1] = potW[6];
            u[2] = potW[5];
            u[3] = potW[0];

            T.mv(u, Tu);

            fluxW[7] = Tu[0];
            potentialDiffW7 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[6];
            u[2] = potNw[5];
            u[3] = potNw[0];

            T.mv(u, Tu);

            fluxNw[7] = Tu[0];
            potentialDiffNw7 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[4];
            u[1] = potW[6];
            u[2] = potW[7];
            u[3] = potW[2];

            T.mv(u, Tu);

            fluxW[7] = Tu[0];
            potentialDiffW7 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[6];
            u[2] = potNw[7];
            u[3] = potNw[2];

            T.mv(u, Tu);

            fluxNw[7] = Tu[0];
            potentialDiffNw7 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[4];
            u[1] = potW[6];
            u[2] = potW[7];
            u[3] = potW[0];

            T.mv(u, Tu);

            fluxW[7] = Tu[0];
            potentialDiffW7 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[6];
            u[2] = potNw[7];
            u[3] = potNw[0];

            T.mv(u, Tu);

            fluxNw[7] = Tu[0];
            potentialDiffNw7 = Tu[0];
        }
        else
        {
            u[0] = potW[4];
            u[1] = potW[6];
            u[2] = potW[5];
            u[3] = potW[2];

            T.mv(u, Tu);

            fluxW[7] = Tu[0];
            potentialDiffW7 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[6];
            u[2] = potNw[5];
            u[3] = potNw[2];

            T.mv(u, Tu);

            fluxNw[7] = Tu[0];
            potentialDiffNw7 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 8)
    {
        // calculate the flux through the subvolumeface 9 (subVolumeFaceIdx = 8)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 4, 0, 6, 2, 5,
                1);

        if (caseL == 1)
        {
            u[0] = potW[4];
            u[1] = potW[0];
            u[2] = potW[6];
            u[3] = potW[5];

            T.mv(u, Tu);

            fluxW[8] = Tu[0];
            potentialDiffW8 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[0];
            u[2] = potNw[6];
            u[3] = potNw[5];

            T.mv(u, Tu);

            fluxNw[8] = Tu[0];
            potentialDiffNw8 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[4];
            u[1] = potW[0];
            u[2] = potW[2];
            u[3] = potW[1];

            T.mv(u, Tu);

            fluxW[8] = Tu[0];
            potentialDiffW8 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[0];
            u[2] = potNw[2];
            u[3] = potNw[1];

            T.mv(u, Tu);

            fluxNw[8] = Tu[0];
            potentialDiffNw8 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[4];
            u[1] = potW[0];
            u[2] = potW[2];
            u[3] = potW[5];

            T.mv(u, Tu);

            fluxW[8] = Tu[0];
            potentialDiffW8 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[0];
            u[2] = potNw[2];
            u[3] = potNw[5];

            T.mv(u, Tu);

            fluxNw[8] = Tu[0];
            potentialDiffNw8 = Tu[0];
        }
        else
        {
            u[0] = potW[4];
            u[1] = potW[0];
            u[2] = potW[6];
            u[3] = potW[1];

            T.mv(u, Tu);

            fluxW[8] = Tu[0];
            potentialDiffW8 = Tu[0];

            u[0] = potNw[4];
            u[1] = potNw[0];
            u[2] = potNw[6];
            u[3] = potNw[1];

            T.mv(u, Tu);

            fluxNw[8] = Tu[0];
            potentialDiffNw8 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 9)
    {
        // calculate the flux through the subvolumeface 10 (subVolumeFaceIdx = 9)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 1, 5, 3, 7, 0,
                4);

        if (caseL == 1)
        {
            u[0] = potW[1];
            u[1] = potW[5];
            u[2] = potW[3];
            u[3] = potW[0];

            T.mv(u, Tu);

            fluxW[9] = Tu[0];
            potentialDiffW9 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[5];
            u[2] = potNw[3];
            u[3] = potNw[0];

            T.mv(u, Tu);

            fluxNw[9] = Tu[0];
            potentialDiffNw9 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[1];
            u[1] = potW[5];
            u[2] = potW[7];
            u[3] = potW[4];

            T.mv(u, Tu);

            fluxW[9] = Tu[0];
            potentialDiffW9 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[5];
            u[2] = potNw[7];
            u[3] = potNw[4];

            T.mv(u, Tu);

            fluxNw[9] = Tu[0];
            potentialDiffNw9 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[1];
            u[1] = potW[5];
            u[2] = potW[7];
            u[3] = potW[0];

            T.mv(u, Tu);

            fluxW[9] = Tu[0];
            potentialDiffW9 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[5];
            u[2] = potNw[7];
            u[3] = potNw[0];

            T.mv(u, Tu);

            fluxNw[9] = Tu[0];
            potentialDiffNw9 = Tu[0];
        }
        else
        {
            u[0] = potW[1];
            u[1] = potW[5];
            u[2] = potW[3];
            u[3] = potW[4];

            T.mv(u, Tu);

            fluxW[9] = Tu[0];
            potentialDiffW9 = Tu[0];

            u[0] = potNw[1];
            u[1] = potNw[5];
            u[2] = potNw[3];
            u[3] = potNw[4];

            T.mv(u, Tu);

            fluxNw[9] = Tu[0];
            potentialDiffNw9 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 10)
    {
        // calculate the flux through the subvolumeface 11 (subVolumeFaceIdx = 10)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 7, 3, 5, 1, 6,
                2);

        if (caseL == 1)
        {
            u[0] = potW[7];
            u[1] = potW[3];
            u[2] = potW[5];
            u[3] = potW[6];

            T.mv(u, Tu);

            fluxW[10] = Tu[0];
            potentialDiffW10 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[3];
            u[2] = potNw[5];
            u[3] = potNw[6];

            T.mv(u, Tu);

            fluxNw[10] = Tu[0];
            potentialDiffNw10 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[7];
            u[1] = potW[3];
            u[2] = potW[1];
            u[3] = potW[2];

            T.mv(u, Tu);

            fluxW[10] = Tu[0];
            potentialDiffW10 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[3];
            u[2] = potNw[1];
            u[3] = potNw[2];

            T.mv(u, Tu);

            fluxNw[10] = Tu[0];
            potentialDiffNw10 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[7];
            u[1] = potW[3];
            u[2] = potW[1];
            u[3] = potW[6];

            T.mv(u, Tu);

            fluxW[10] = Tu[0];
            potentialDiffW10 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[3];
            u[2] = potNw[1];
            u[3] = potNw[6];

            T.mv(u, Tu);

            fluxNw[10] = Tu[0];
            potentialDiffNw10 = Tu[0];
        }
        else
        {
            u[0] = potW[7];
            u[1] = potW[3];
            u[2] = potW[5];
            u[3] = potW[2];

            T.mv(u, Tu);

            fluxW[10] = Tu[0];
            potentialDiffW10 = Tu[0];

            u[0] = potNw[7];
            u[1] = potNw[3];
            u[2] = potNw[5];
            u[3] = potNw[2];

            T.mv(u, Tu);

            fluxNw[10] = Tu[0];
            potentialDiffNw10 = Tu[0];
        }
    }

    if (fIdx < 0 || fIdx == 11)
    {
        // calculate the flux through the subvolumeface 12 (subVolumeFaceIdx = 11)
        int caseL = transmissibilityCalculator.transmissibility(T, interactionVolume, lambda, 2, 6, 0, 4, 3,
                7);

        if (caseL == 1)
        {
            u[0] = potW[2];
            u[1] = potW[6];
            u[2] = potW[0];
            u[3] = potW[3];

            T.mv(u, Tu);

            fluxW[11] = Tu[0];
            potentialDiffW11 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[6];
            u[2] = potNw[0];
            u[3] = potNw[3];

            T.mv(u, Tu);

            fluxNw[11] = Tu[0];
            potentialDiffNw11 = Tu[0];
        }
        else if (caseL == 2)
        {
            u[0] = potW[2];
            u[1] = potW[6];
            u[2] = potW[4];
            u[3] = potW[7];

            T.mv(u, Tu);

            fluxW[11] = Tu[0];
            potentialDiffW11 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[6];
            u[2] = potNw[4];
            u[3] = potNw[7];

            T.mv(u, Tu);

            fluxNw[11] = Tu[0];
            potentialDiffNw11 = Tu[0];
        }
        else if (caseL == 3)
        {
            u[0] = potW[2];
            u[1] = potW[6];
            u[2] = potW[4];
            u[3] = potW[3];

            T.mv(u, Tu);

            fluxW[11] = Tu[0];
            potentialDiffW11 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[6];
            u[2] = potNw[4];
            u[3] = potNw[3];

            T.mv(u, Tu);

            fluxNw[11] = Tu[0];
            potentialDiffNw11 = Tu[0];
        }
        else
        {
            u[0] = potW[2];
            u[1] = potW[6];
            u[2] = potW[0];
            u[3] = potW[7];

            T.mv(u, Tu);

            fluxW[11] = Tu[0];
            potentialDiffW11 = Tu[0];

            u[0] = potNw[2];
            u[1] = potNw[6];
            u[2] = potNw[0];
            u[3] = potNw[7];

            T.mv(u, Tu);

            fluxNw[11] = Tu[0];
            potentialDiffNw11 = Tu[0];
        }
    }

    //store potentials for further calculations (saturation, ...)
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialDiffW0);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialDiffNw0);
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 1), -potentialDiffW3);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 1), -potentialDiffNw3);
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 2), -potentialDiffW8);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 2), -potentialDiffNw8);

    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 0), potentialDiffW1);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 0), potentialDiffNw1);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialDiffW0);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialDiffNw0);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 2), potentialDiffW9);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 2), potentialDiffNw9);

    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 0), potentialDiffW3);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 0), potentialDiffNw3);
    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 1), -potentialDiffW2);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 1), -potentialDiffNw2);
    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 2), potentialDiffW11);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 2), potentialDiffNw11);

    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 0), potentialDiffW2);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 0), potentialDiffNw2);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialDiffW1);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialDiffNw1);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 2), -potentialDiffW10);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 2), -potentialDiffNw10);

    cellData5.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(4, 0), potentialDiffW8);
    cellData5.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(4, 0), potentialDiffNw8);
    cellData5.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(4, 1), -potentialDiffW4);
    cellData5.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(4, 1), -potentialDiffNw4);
    cellData5.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(4, 2), potentialDiffW7);
    cellData5.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(4, 2), potentialDiffNw7);

    cellData6.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(5, 0), -potentialDiffW9);
    cellData6.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(5, 0), -potentialDiffNw9);
    cellData6.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(5, 1), -potentialDiffW5);
    cellData6.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(5, 1), -potentialDiffNw5);
    cellData6.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(5, 2), potentialDiffW4);
    cellData6.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(5, 2), potentialDiffNw4);

    cellData7.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(6, 0), -potentialDiffW11);
    cellData7.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(6, 0), -potentialDiffNw11);
    cellData7.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(6, 1), -potentialDiffW7);
    cellData7.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(6, 1), -potentialDiffNw7);
    cellData7.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(6, 2), potentialDiffW6);
    cellData7.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(6, 2), potentialDiffNw6);

    cellData8.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(7, 0), potentialDiffW10);
    cellData8.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(7, 0), potentialDiffNw10);
    cellData8.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(7, 1), -potentialDiffW6);
    cellData8.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(7, 1), -potentialDiffNw6);
    cellData8.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(7, 2), potentialDiffW5);
    cellData8.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(7, 2), potentialDiffNw5);

    // compute mobilities of subvolumeface 1 (subVolumeFaceIdx = 0)
    Dune::FieldVector<Scalar, numPhases> lambda0Upw(0.0);
    lambda0Upw[wPhaseIdx] = (potentialDiffW0 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
    lambda0Upw[nPhaseIdx] = (potentialDiffNw0 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

    // compute mobilities of subvolumeface 2 (subVolumeFaceIdx = 1)
    Dune::FieldVector<Scalar, numPhases> lambda1Upw(0.0);
    lambda1Upw[wPhaseIdx] = (potentialDiffW1 >= 0) ? lambda2[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda1Upw[nPhaseIdx] = (potentialDiffNw1 >= 0) ? lambda2[nPhaseIdx] : lambda4[nPhaseIdx];

    // compute mobilities of subvolumeface 3 (subVolumeFaceIdx = 2)
    Dune::FieldVector<Scalar, numPhases> lambda2Upw(0.0);
    lambda2Upw[wPhaseIdx] = (potentialDiffW2 >= 0) ? lambda4[wPhaseIdx] : lambda3[wPhaseIdx];
    lambda2Upw[nPhaseIdx] = (potentialDiffNw2 >= 0) ? lambda4[nPhaseIdx] : lambda3[nPhaseIdx];

    // compute mobilities of subvolumeface 4 (subVolumeFaceIdx = 3)
    Dune::FieldVector<Scalar, numPhases> lambda3Upw(0.0);
    lambda3Upw[wPhaseIdx] = (potentialDiffW3 >= 0) ? lambda3[wPhaseIdx] : lambda1[wPhaseIdx];
    lambda3Upw[nPhaseIdx] = (potentialDiffNw3 >= 0) ? lambda3[nPhaseIdx] : lambda1[nPhaseIdx];

    // compute mobilities of subvolumeface 5 (subVolumeFaceIdx = 4)
    Dune::FieldVector<Scalar, numPhases> lambda4Upw(0.0);
    lambda4Upw[wPhaseIdx] = (potentialDiffW4 >= 0) ? lambda6[wPhaseIdx] : lambda5[wPhaseIdx];
    lambda4Upw[nPhaseIdx] = (potentialDiffNw4 >= 0) ? lambda6[nPhaseIdx] : lambda5[nPhaseIdx];

    // compute mobilities of subvolumeface 6 (subVolumeFaceIdx = 5)
    Dune::FieldVector<Scalar, numPhases> lambda5Upw(0.0);
    lambda5Upw[wPhaseIdx] = (potentialDiffW5 >= 0) ? lambda8[wPhaseIdx] : lambda6[wPhaseIdx];
    lambda5Upw[nPhaseIdx] = (potentialDiffNw5 >= 0) ? lambda8[nPhaseIdx] : lambda6[nPhaseIdx];

    // compute mobilities of subvolumeface 7 (subVolumeFaceIdx = 6)
    Dune::FieldVector<Scalar, numPhases> lambda6Upw(0.0);
    lambda6Upw[wPhaseIdx] = (potentialDiffW6 >= 0) ? lambda7[wPhaseIdx] : lambda8[wPhaseIdx];
    lambda6Upw[nPhaseIdx] = (potentialDiffNw6 >= 0) ? lambda7[nPhaseIdx] : lambda8[nPhaseIdx];

    // compute mobilities of subvolumeface 8 (subVolumeFaceIdx = 7)
    Dune::FieldVector<Scalar, numPhases> lambda7Upw(0.0);
    lambda7Upw[wPhaseIdx] = (potentialDiffW7 >= 0) ? lambda5[wPhaseIdx] : lambda7[wPhaseIdx];
    lambda7Upw[nPhaseIdx] = (potentialDiffNw7 >= 0) ? lambda5[nPhaseIdx] : lambda7[nPhaseIdx];

    // compute mobilities of subvolumeface 9 (subVolumeFaceIdx = 8)
    Dune::FieldVector<Scalar, numPhases> lambda8Upw(0.0);
    lambda8Upw[wPhaseIdx] = (potentialDiffW8 >= 0) ? lambda5[wPhaseIdx] : lambda1[wPhaseIdx];
    lambda8Upw[nPhaseIdx] = (potentialDiffNw8 >= 0) ? lambda5[nPhaseIdx] : lambda1[nPhaseIdx];

    // compute mobilities of subvolumeface 10 (subVolumeFaceIdx = 9)
    Dune::FieldVector<Scalar, numPhases> lambda9Upw(0.0);
    lambda9Upw[wPhaseIdx] = (potentialDiffW9 >= 0) ? lambda2[wPhaseIdx] : lambda6[wPhaseIdx];
    lambda9Upw[nPhaseIdx] = (potentialDiffNw9 >= 0) ? lambda2[nPhaseIdx] : lambda6[nPhaseIdx];

    // compute mobilities of subvolumeface 11 (subVolumeFaceIdx = 10)
    Dune::FieldVector<Scalar, numPhases> lambda10Upw(0.0);
    lambda10Upw[wPhaseIdx] = (potentialDiffW10 >= 0) ? lambda8[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda10Upw[nPhaseIdx] = (potentialDiffNw10 >= 0) ? lambda8[nPhaseIdx] : lambda4[nPhaseIdx];

    // compute mobilities of subvolumeface 12 (subVolumeFaceIdx = 11)
    Dune::FieldVector<Scalar, numPhases> lambda11Upw(0.0);
    lambda11Upw[wPhaseIdx] = (potentialDiffW11 >= 0) ? lambda3[wPhaseIdx] : lambda7[wPhaseIdx];
    lambda11Upw[nPhaseIdx] = (potentialDiffNw11 >= 0) ? lambda3[nPhaseIdx] : lambda7[nPhaseIdx];

    for (int i = 0; i < numPhases; i++)
    {
        // evaluate parts of velocity --> always take the normal for which the flux is calculated!
        DimVector vel12 = interactionVolume.getNormal(0, 0);
        DimVector vel21 = interactionVolume.getNormal(1, 1);
        DimVector vel24 = interactionVolume.getNormal(1, 0);
        DimVector vel42 = interactionVolume.getNormal(3, 1);
        DimVector vel43 = interactionVolume.getNormal(3, 0);
        DimVector vel34 = interactionVolume.getNormal(2, 1);
        DimVector vel31 = interactionVolume.getNormal(2, 0);
        DimVector vel13 = interactionVolume.getNormal(0, 1);
        DimVector vel65 = interactionVolume.getNormal(5, 2);
        DimVector vel56 = interactionVolume.getNormal(4, 1);
        DimVector vel86 = interactionVolume.getNormal(7, 2);
        DimVector vel68 = interactionVolume.getNormal(5, 1);
        DimVector vel78 = interactionVolume.getNormal(6, 2);
        DimVector vel87 = interactionVolume.getNormal(7, 1);
        DimVector vel57 = interactionVolume.getNormal(4, 2);
        DimVector vel75 = interactionVolume.getNormal(6, 1);
        DimVector vel51 = interactionVolume.getNormal(4, 0);
        DimVector vel15 = interactionVolume.getNormal(0, 2);
        DimVector vel26 = interactionVolume.getNormal(1, 2);
        DimVector vel62 = interactionVolume.getNormal(5, 0);
        DimVector vel84 = interactionVolume.getNormal(7, 0);
        DimVector vel48 = interactionVolume.getNormal(3, 2);
        DimVector vel37 = interactionVolume.getNormal(2, 2);
        DimVector vel73 = interactionVolume.getNormal(6, 0);
        vel21 *= -1;
        vel42 *= -1;
        vel34 *= -1;
        vel13 *= -1;
        vel56 *= -1;
        vel68 *= -1;
        vel87 *= -1;
        vel75 *= -1;
        vel15 *= -1;
        vel62 *= -1;
        vel48 *= -1;
        vel73 *= -1;

        Dune::FieldVector<Scalar, 12> flux(0);
        switch (i)
        {
        case wPhaseIdx:
        {
            flux = fluxW;
            break;
        }
        case nPhaseIdx:
        {
            flux = fluxNw;
            break;
        }
        }

        vel12 *= flux[0] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal1, 0, 0));
        vel21 *= flux[0] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal2, 1, 1));
        vel24 *= flux[1] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal2, 1, 0));
        vel42 *= flux[1] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal4, 3, 1));
        vel43 *= flux[2] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal4, 3, 0));
        vel34 *= flux[2] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal3, 2, 1));
        vel31 *= flux[3] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal3, 2, 0));
        vel13 *= flux[3] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal1, 0, 1));
        vel65 *= flux[4] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal6, 5, 2));
        vel56 *= flux[4] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal5, 4, 1));
        vel86 *= flux[5] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal8, 7, 2));
        vel68 *= flux[5] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal6, 5, 1));
        vel78 *= flux[6] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal7, 6, 2));
        vel87 *= flux[6] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal8, 7, 1));
        vel57 *= flux[7] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal5, 4, 2));
        vel75 *= flux[7] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal7, 6, 1));
        vel51 *= flux[8] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal5, 4, 0));
        vel15 *= flux[8] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal1, 0, 2));
        vel26 *= flux[9] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal2, 1, 2));
        vel62 *= flux[9] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal6, 5, 0));
        vel84 *= flux[10] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal8, 7, 0));
        vel48 *= flux[10] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal4, 3, 2));
        vel37 *= flux[11] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal3, 2, 2));
        vel73 *= flux[11] / (interactionVolumes.getRealFluxFaceArea(interactionVolume, eIdxGlobal7, 6, 0));


        Scalar lambdaT0 = lambda0Upw[wPhaseIdx] + lambda0Upw[nPhaseIdx];
        Scalar lambdaT1 = lambda1Upw[wPhaseIdx] + lambda1Upw[nPhaseIdx];
        Scalar lambdaT2 = lambda2Upw[wPhaseIdx] + lambda2Upw[nPhaseIdx];
        Scalar lambdaT3 = lambda3Upw[wPhaseIdx] + lambda3Upw[nPhaseIdx];
        Scalar lambdaT4 = lambda4Upw[wPhaseIdx] + lambda4Upw[nPhaseIdx];
        Scalar lambdaT5 = lambda5Upw[wPhaseIdx] + lambda5Upw[nPhaseIdx];
        Scalar lambdaT6 = lambda6Upw[wPhaseIdx] + lambda6Upw[nPhaseIdx];
        Scalar lambdaT7 = lambda7Upw[wPhaseIdx] + lambda7Upw[nPhaseIdx];
        Scalar lambdaT8 = lambda8Upw[wPhaseIdx] + lambda8Upw[nPhaseIdx];
        Scalar lambdaT9 = lambda9Upw[wPhaseIdx] + lambda9Upw[nPhaseIdx];
        Scalar lambdaT10 = lambda10Upw[wPhaseIdx] + lambda10Upw[nPhaseIdx];
        Scalar lambdaT11 = lambda11Upw[wPhaseIdx] + lambda11Upw[nPhaseIdx];
        Scalar fracFlow0 = (lambdaT0 > threshold_) ? lambda0Upw[i] / (lambdaT0) : 0.0;
        Scalar fracFlow1 = (lambdaT1 > threshold_) ? lambda1Upw[i] / (lambdaT1) : 0.0;
        Scalar fracFlow2 = (lambdaT2 > threshold_) ? lambda2Upw[i] / (lambdaT2) : 0.0;
        Scalar fracFlow3 = (lambdaT3 > threshold_) ? lambda3Upw[i] / (lambdaT3) : 0.0;
        Scalar fracFlow4 = (lambdaT4 > threshold_) ? lambda4Upw[i] / (lambdaT4) : 0.0;
        Scalar fracFlow5 = (lambdaT5 > threshold_) ? lambda5Upw[i] / (lambdaT5) : 0.0;
        Scalar fracFlow6 = (lambdaT6 > threshold_) ? lambda6Upw[i] / (lambdaT6) : 0.0;
        Scalar fracFlow7 = (lambdaT7 > threshold_) ? lambda7Upw[i] / (lambdaT7) : 0.0;
        Scalar fracFlow8 = (lambdaT8 > threshold_) ? lambda8Upw[i] / (lambdaT8) : 0.0;
        Scalar fracFlow9 = (lambdaT9 > threshold_) ? lambda9Upw[i] / (lambdaT9) : 0.0;
        Scalar fracFlow10 = (lambdaT10 > threshold_) ? lambda10Upw[i] / (lambdaT10) : 0.0;
        Scalar fracFlow11 = (lambdaT11 > threshold_) ? lambda11Upw[i] / (lambdaT11) : 0.0;

        vel12 *= fracFlow0;
        vel21 *= fracFlow0;
        vel24 *= fracFlow1;
        vel42 *= fracFlow1;
        vel43 *= fracFlow2;
        vel34 *= fracFlow2;
        vel31 *= fracFlow3;
        vel13 *= fracFlow3;
        vel65 *= fracFlow4;
        vel56 *= fracFlow4;
        vel86 *= fracFlow5;
        vel68 *= fracFlow5;
        vel78 *= fracFlow6;
        vel87 *= fracFlow6;
        vel57 *= fracFlow7;
        vel75 *= fracFlow7;
        vel51 *= fracFlow8;
        vel15 *= fracFlow8;
        vel26 *= fracFlow9;
        vel62 *= fracFlow9;
        vel84 *= fracFlow10;
        vel48 *= fracFlow10;
        vel37 *= fracFlow11;
        vel73 *= fracFlow11;

        //store velocities
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 0), vel12);
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 1), vel13);
        cellData1.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(0, 2), vel15);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 0), vel24);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 1), vel21);
        cellData2.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(1, 2), vel26);
        cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 0), vel31);
        cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 1), vel34);
        cellData3.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(2, 2), vel37);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 0), vel43);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 1), vel42);
        cellData4.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(3, 2), vel48);
        cellData5.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(4, 0), vel51);
        cellData5.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(4, 1), vel56);
        cellData5.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(4, 2), vel57);
        cellData6.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(5, 0), vel62);
        cellData6.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(5, 1), vel68);
        cellData6.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(5, 2), vel65);
        cellData7.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(6, 0), vel73);
        cellData7.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(6, 1), vel75);
        cellData7.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(6, 2), vel78);
        cellData8.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(7, 0), vel84);
        cellData8.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(7, 1), vel87);
        cellData8.fluxData().addVelocity(i, interactionVolume.getIndexOnElement(7, 2), vel86);
    }
    //set velocity marker
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 0));
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 1));
    cellData1.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(0, 2));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 0));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 1));
    cellData2.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(1, 2));
    cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 0));
    cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 1));
    cellData3.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(2, 2));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 0));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 1));
    cellData4.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(3, 2));
    cellData5.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(4, 0));
    cellData5.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(4, 1));
    cellData5.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(4, 2));
    cellData6.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(5, 0));
    cellData6.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(5, 1));
    cellData6.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(5, 2));
    cellData7.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(6, 0));
    cellData7.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(6, 1));
    cellData7.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(6, 2));
    cellData8.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(7, 0));
    cellData8.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(7, 1));
    cellData8.fluxData().setVelocityMarker(interactionVolume.getIndexOnElement(7, 2));
        }

/*!
 * \brief Calculates the velocity at a boundary flux faces.
 *
 *  Calculates the velocity at a boundary flux face and adds it to the face velocity vector in the <tt>CellData</tt> object.
 *
 * \param interactionVolume An <tt>InteractionVolume</tt> object including the information for calculating the MPFA transmissibilities
 * \param cellData  <tt>CellData</tt> object of an IMPES model for the sub-volume including the boundary face
 * \param elemIdx local sub-volume index
 */
template<class TypeTag>
void FvMpfaL3dVelocity2p<TypeTag>::calculateBoundaryInteractionVolumeVelocity(InteractionVolume& interactionVolume, CellData& cellData, int elemIdx)
{
    auto element = interactionVolume.getSubVolumeElement(elemIdx);

    // get global coordinate of cell centers
    const GlobalPosition& globalPos = element.geometry().center();

    // permeability vector at boundary
    DimMatrix permeability(problem_.spatialParams().intrinsicPermeability(element));

    //get mobilities of the phases
    Dune::FieldVector<Scalar, numPhases> lambda(cellData.mobility(wPhaseIdx));
    lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

    Scalar potW = cellData.potential(wPhaseIdx);
    Scalar potNw = cellData.potential(nPhaseIdx);

    for (int fIdx = 0; fIdx < dim; fIdx++)
    {
        int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);

        if (interactionVolume.isBoundaryFace(intVolFaceIdx))
        {
            if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(pressureEqIdx))
            {
                int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, fIdx);

                const GlobalPosition& globalPosFace = interactionVolume.getFacePosition(elemIdx, fIdx);

                DimVector distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();
                DimVector& normal = interactionVolume.getNormal(elemIdx, fIdx);

                // get pc and lambda at the boundary
                Scalar satWBound = cellData.saturation(wPhaseIdx);
                //check boundary sat at face 1
                if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(satEqIdx))
                {
                    Scalar satBound = interactionVolume.getDirichletValues(intVolFaceIdx)[saturationIdx];
                    switch (saturationType_)
                    {
                    case sw:
                    {
                        satWBound = satBound;
                        break;
                    }
                    case sn:
                    {
                        satWBound = 1 - satBound;
                        break;
                    }
                    }

                }

                // old material law interface is deprecated: Replace this by
                // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
                // after the release of 3.3, when the deprecated interface is no longer supported
                const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

                Scalar pcBound = fluidMatrixInteraction.pc(satWBound);

                Scalar gravityDiffBound = (problem_.bBoxMax() - globalPosFace) * gravity_
                        * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                pcBound += gravityDiffBound;

                Dune::FieldVector<Scalar, numPhases> lambdaBound(fluidMatrixInteraction.krw(satWBound));
                lambdaBound[nPhaseIdx] = fluidMatrixInteraction.krn(satWBound);
                lambdaBound[wPhaseIdx] /= viscosity_[wPhaseIdx];
                lambdaBound[nPhaseIdx] /= viscosity_[nPhaseIdx];

                Scalar gdeltaZ = (problem_.bBoxMax()-globalPosFace) * gravity_;
                Scalar potentialBoundW = interactionVolume.getDirichletValues(intVolFaceIdx)[pressureIdx] + density_[wPhaseIdx]*gdeltaZ;
                Scalar potentialBoundNw = potentialBoundW;

                //calculate potential gradients
                switch (pressureType_)
                {
                case pw:
                {
                    potentialBoundNw += pcBound;
                    break;
                }
                case pn:
                {
                    //calculate potential gradients
                    potentialBoundW -= pcBound;
                    break;
                }
                }

                Scalar potentialDiffW = (potW - potentialBoundW) / dist;
                Scalar potentialDiffNw = (potNw - potentialBoundNw) / dist;

                //store potentials for further calculations (saturation, ...)
                cellData.fluxData().addUpwindPotential(wPhaseIdx, boundaryFaceIdx, potentialDiffW);
                cellData.fluxData().addUpwindPotential(nPhaseIdx, boundaryFaceIdx, potentialDiffNw);

                //calculated phase velocities from advective velocities -> capillary pressure velocity already added in pressure part!
                DimVector velocityW(0);
                DimVector velocityNw(0);

                // calculate capillary pressure gradient
                DimVector potentialGradient = normal;
                potentialGradient *= (potW - potentialBoundW) / dist;
                permeability.mv(potentialGradient, velocityW);

                potentialGradient = normal;
                potentialGradient *= (potNw - potentialBoundNw) / dist;
                permeability.mv(potentialGradient, velocityNw);

                velocityW *= (potentialDiffW >= 0.) ? lambda[wPhaseIdx] : lambdaBound[wPhaseIdx];
                velocityNw *= (potentialDiffNw >= 0.) ? lambda[nPhaseIdx] : lambdaBound[nPhaseIdx];

                //velocity is calculated from two vertices of one intersection!
                velocityW *= 0.25;
                velocityNw *= 0.25;

                //store velocities
                cellData.fluxData().addVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                cellData.fluxData().addVelocity(nPhaseIdx, boundaryFaceIdx, velocityNw);
                cellData.fluxData().setVelocityMarker(boundaryFaceIdx);
            }
            else if (interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressureEqIdx))
            {
                int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, fIdx);

                DimVector& normal = interactionVolume.getNormal(elemIdx, fIdx);

                // get neumann boundary value
                PrimaryVariables boundValues(interactionVolume.getNeumannValues(intVolFaceIdx));

                boundValues[wPhaseIdx] /= density_[wPhaseIdx];
                boundValues[nPhaseIdx] /= density_[nPhaseIdx];

                DimVector velocityW(normal);
                DimVector velocityNw(normal);

                velocityW *= boundValues[wPhaseIdx] / (4.0*interactionVolume.getFaceArea(elemIdx, fIdx));
                velocityNw *= boundValues[nPhaseIdx]
                                          / (4.0*interactionVolume.getFaceArea(elemIdx, fIdx));

                //store potentials for further calculations (saturation, ...)
                cellData.fluxData().addUpwindPotential(wPhaseIdx, boundaryFaceIdx, boundValues[wPhaseIdx]);
                cellData.fluxData().addUpwindPotential(nPhaseIdx, boundaryFaceIdx, boundValues[nPhaseIdx]);

                //store velocities
                cellData.fluxData().addVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                cellData.fluxData().addVelocity(nPhaseIdx, boundaryFaceIdx, velocityNw);
                cellData.fluxData().setVelocityMarker(boundaryFaceIdx);
            }
            else
            {
                std::cout << "interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressureEqIdx)"
                        << interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressureEqIdx) << "\n";
                DUNE_THROW(Dune::NotImplemented,
                        "No valid boundary condition type defined for pressure equation!");
            }
        }
    }
}

} // end namespace Dumux
#endif
