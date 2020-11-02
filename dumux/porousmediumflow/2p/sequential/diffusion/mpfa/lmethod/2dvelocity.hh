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
 * \brief  Velocity calculation using a 2-d MPFA L-method.
 */
#ifndef DUMUX_FVMPFAL2DVELOCITY2P_HH
#define DUMUX_FVMPFAL2DVELOCITY2P_HH

#include <dune/grid/common/gridenums.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/linteractionvolume.hh>
#include "2dtransmissibilitycalculator.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Class for calculating 2-d velocities from cell-wise constant pressure values.
 *
 * Calculates phase velocities or total velocity from a known pressure field applying a
 * finite volume discretization and a MPFA L-method.
 * At Dirichlet boundaries a two-point flux approximation is used.
 * The pressure has to be given as piecewise constant cell values.
 * The velocities are calculated as
 *
 *\f[ \boldsymbol v_\alpha = - \lambda_\alpha \boldsymbol K \textbf{grad}\, \Phi_\alpha, \f]
 * and,
 * \f[ \boldsymbol v_t = \boldsymbol v_w + \boldsymbol v_n,\f]
 *
 * where \f$ \Phi_\alpha \f$ denotes the potential of phase \f$ \alpha \f$,
 * \f$ \boldsymbol K \f$ the intrinsic permeability,
 * and \f$ \lambda_\alpha \f$ a phase mobility.
 *
 * Remark1: only for quadrilaterals!
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag> class FvMpfaL2dVelocity2p
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
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Grid = typename GridView::Grid;
    using IndexSet = typename GridView::IndexSet;

    using Geometry = typename Element::Geometry;
    using JacobianTransposed = typename Geometry::JacobianTransposed ;

    using GridTypeIndices = GetPropType<TypeTag, Properties::GridTypeIndices>;

    using InteractionVolume = FVMPFALInteractionVolume<TypeTag>;
    using TransmissibilityCalculator = FvMpfaL2dTransmissibilityCalculator<TypeTag>;
    using InnerBoundaryVolumeFaces = std::vector<Dune::FieldVector<bool, 2*dim> >;

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

    using LocalPosition = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    /*!
     * \brief Constructs a FvMpfaL2dVelocity2p object
     *
     * \param problem A problem class object
     */
    FvMpfaL2dVelocity2p(Problem& problem) :
        problem_(problem), transmissibilityCalculator_(problem), gravity_(problem.gravity())
    {
        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        vtkOutputLevel_ = getParam<int>("Vtk.OutputLevel");
    }

    //! Calculate velocities for flux faces of an interaction volume
    void calculateInnerInteractionVolumeVelocity(InteractionVolume& interactionVolume,
                            CellData& cellData1, CellData& cellData2, CellData& cellData3, CellData& cellData4,
                            InnerBoundaryVolumeFaces& innerBoundaryVolumeFaces);
    void calculateBoundaryInteractionVolumeVelocity(InteractionVolume& interactionVolume,
                            CellData& cellData, int elemIdx);

    //! Initializes the velocity model
    void initialize()
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
            Dune::BlockVector < DimVector > &velocityWetting
                    = *(writer.template allocateManagedBuffer<Scalar,dim>(problem_.gridView().size(0)));
            Dune::BlockVector < DimVector > &velocityNonwetting
                    = *(writer.template allocateManagedBuffer<Scalar,dim>(problem_.gridView().size(0)));

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
protected:
    TransmissibilityCalculator transmissibilityCalculator_;

    const GravityVector& gravity_; // vector including the gravity constant

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
 *  Calculates the velocities at the flux faces of an interation volume around a vertex which is not
 *  a boundary vertex and adds them to the face velocity vectors in the <tt>CellData</tt> objects.
 *
 * \param interactionVolume An <tt>InteractionVolume</tt> object including the information for calculating the MPFA transmissibilities
 * \param cellData1  <tt>CellData</tt> object of an IMPES model for sub-volume 1
 * \param cellData2  <tt>CellData</tt> object of an IMPES model for sub-volume 2
 * \param cellData3  <tt>CellData</tt> object of an IMPES model for sub-volume 3
 * \param cellData4  <tt>CellData</tt> object of an IMPES model for sub-volume 4
 * \param innerBoundaryVolumeFaces container including information about faces intersecting a boundary
 */
template<class TypeTag>
void FvMpfaL2dVelocity2p<TypeTag>::calculateInnerInteractionVolumeVelocity(InteractionVolume& interactionVolume,
                                                                           CellData& cellData1, CellData& cellData2,
                                                                           CellData& cellData3, CellData& cellData4,
                                                                           InnerBoundaryVolumeFaces& innerBoundaryVolumeFaces)
{
    auto element1 = interactionVolume.getSubVolumeElement(0);
    auto element2 = interactionVolume.getSubVolumeElement(1);
    auto element3 = interactionVolume.getSubVolumeElement(2);
    auto element4 = interactionVolume.getSubVolumeElement(3);

    int level1 = element1.level();
    int level2 = element2.level();
    int level3 = element3.level();
    int level4 = element4.level();

    // cell index
    int eIdxGlobal1 = problem_.variables().index(element1);
    int eIdxGlobal2 = problem_.variables().index(element2);
    int eIdxGlobal3 = problem_.variables().index(element3);
    int eIdxGlobal4 = problem_.variables().index(element4);

    // get pressure values
    Dune::FieldVector < Scalar, 2 * dim > potW(0);
    Dune::FieldVector < Scalar, 2 * dim > potNw(0);

    potW[0] = cellData1.potential(wPhaseIdx);
    potW[1] = cellData2.potential(wPhaseIdx);
    potW[2] = cellData3.potential(wPhaseIdx);
    potW[3] = cellData4.potential(wPhaseIdx);

    potNw[0] = cellData1.potential(nPhaseIdx);
    potNw[1] = cellData2.potential(nPhaseIdx);
    potNw[2] = cellData3.potential(nPhaseIdx);
    potNw[3] = cellData4.potential(nPhaseIdx);

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda1(cellData1.mobility(wPhaseIdx));
    lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda2(cellData2.mobility(wPhaseIdx));
    lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda3(cellData3.mobility(wPhaseIdx));
    lambda3[nPhaseIdx] = cellData3.mobility(nPhaseIdx);

    //compute total mobility of cell 1
    Scalar lambdaTotal3 = lambda3[wPhaseIdx] + lambda3[nPhaseIdx];

    //get mobilities of the phases
    Dune::FieldVector < Scalar, numPhases > lambda4(cellData4.mobility(wPhaseIdx));
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

    Scalar potentialDiffW12 = 0;
    Scalar potentialDiffW14 = 0;
    Scalar potentialDiffW32 = 0;
    Scalar potentialDiffW34 = 0;

    Scalar potentialDiffNw12 = 0;
    Scalar potentialDiffNw14 = 0;
    Scalar potentialDiffNw32 = 0;
    Scalar potentialDiffNw34 = 0;

    //flux vector
    Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
    Dune::FieldVector < Scalar, 2 * dim > fluxNw(0);

    Dune::FieldMatrix < Scalar, dim, 2 * dim - dim + 1 > T(0);
    DimVector Tu(0);
    Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);

    int lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 0, 1, 2, 3);

    if (lType == TransmissibilityCalculator::rightTriangle)
    {
        u[0] = potW[1];
        u[1] = potW[2];
        u[2] = potW[0];

        T.mv(u, Tu);

        fluxW[0] = Tu[1];
        potentialDiffW12 = Tu[1];

        u[0] = potNw[1];
        u[1] = potNw[2];
        u[2] = potNw[0];

        T.mv(u, Tu);

        fluxNw[0] = Tu[1];
        potentialDiffNw12 = Tu[1];
    }
    else
    {
        u[0] = potW[0];
        u[1] = potW[3];
        u[2] = potW[1];

        T.mv(u, Tu);

        fluxW[0] = Tu[1];
        potentialDiffW12 = Tu[1];

        u[0] = potNw[0];
        u[1] = potNw[3];
        u[2] = potNw[1];

        T.mv(u, Tu);

        fluxNw[0] = Tu[1];
        potentialDiffNw12 = Tu[1];
    }

    lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 1, 2, 3, 0);

    if (lType == TransmissibilityCalculator::rightTriangle)
    {
        u[0] = potW[2];
        u[1] = potW[3];
        u[2] = potW[1];

        T.mv(u, Tu);

        fluxW[1] = Tu[1];
        potentialDiffW32 = -Tu[1];

        u[0] = potNw[2];
        u[1] = potNw[3];
        u[2] = potNw[1];

        T.mv(u, Tu);

        fluxNw[1] = Tu[1];
        potentialDiffNw32 = -Tu[1];
    }
    else
    {
        u[0] = potW[1];
        u[1] = potW[0];
        u[2] = potW[2];

        T.mv(u, Tu);

        fluxW[1] = Tu[1];
        potentialDiffW32 = -Tu[1];

        u[0] = potNw[1];
        u[1] = potNw[0];
        u[2] = potNw[2];

        T.mv(u, Tu);

        fluxNw[1] = Tu[1];
        potentialDiffNw32 = -Tu[1];
    }

    lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 2, 3, 0, 1);

    if (lType == TransmissibilityCalculator::rightTriangle)
    {
        u[0] = potW[3];
        u[1] = potW[0];
        u[2] = potW[2];

        T.mv(u, Tu);

        fluxW[2] = Tu[1];
        potentialDiffW34 = Tu[1];

        u[0] = potNw[3];
        u[1] = potNw[0];
        u[2] = potNw[2];

        T.mv(u, Tu);

        fluxNw[2] = Tu[1];
        potentialDiffNw34 = Tu[1];
    }
    else
    {
        u[0] = potW[2];
        u[1] = potW[1];
        u[2] = potW[3];

        T.mv(u, Tu);

        fluxW[2] = Tu[1];
        potentialDiffW34 = Tu[1];

        u[0] = potNw[2];
        u[1] = potNw[1];
        u[2] = potNw[3];

        T.mv(u, Tu);

        fluxNw[2] = Tu[1];
        potentialDiffNw34 = Tu[1];
    }

    lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 3, 0, 1, 2);

    if (lType == TransmissibilityCalculator::rightTriangle)
    {
        u[0] = potW[0];
        u[1] = potW[1];
        u[2] = potW[3];

        T.mv(u, Tu);

        fluxW[3] = Tu[1];
        potentialDiffW14 = -Tu[1];

        u[0] = potNw[0];
        u[1] = potNw[1];
        u[2] = potNw[3];

        T.mv(u, Tu);

        fluxNw[3] = Tu[1];
        potentialDiffNw14 = -Tu[1];
    }
    else
    {
        u[0] = potW[3];
        u[1] = potW[2];
        u[2] = potW[0];

        T.mv(u, Tu);

        fluxW[3] = Tu[1];
        potentialDiffW14 = -Tu[1];

        u[0] = potNw[3];
        u[1] = potNw[2];
        u[2] = potNw[0];

        T.mv(u, Tu);

        fluxNw[3] = Tu[1];
        potentialDiffNw14 = -Tu[1];
    }

    //store potentials for further calculations (saturation, ...)
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialDiffW12);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 0), potentialDiffNw12);
    cellData1.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(0, 1), potentialDiffW14);
    cellData1.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(0, 1), potentialDiffNw14);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 0), -potentialDiffW32);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 0), -potentialDiffNw32);
    cellData2.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialDiffW12);
    cellData2.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(1, 1), -potentialDiffNw12);
    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 0), potentialDiffW34);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 0), potentialDiffNw34);
    cellData3.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(2, 1), potentialDiffW32);
    cellData3.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(2, 1), potentialDiffNw32);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -potentialDiffW14);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 0), -potentialDiffNw14);
    cellData4.fluxData().addUpwindPotential(wPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialDiffW34);
    cellData4.fluxData().addUpwindPotential(nPhaseIdx, interactionVolume.getIndexOnElement(3, 1), -potentialDiffNw34);

    //compute mobilities of face 1
    Dune::FieldVector < Scalar, numPhases > lambda12Upw(0.0);
    lambda12Upw[wPhaseIdx] = (potentialDiffW12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
    lambda12Upw[nPhaseIdx] = (potentialDiffNw12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

    //compute mobilities of face 4
    Dune::FieldVector < Scalar, numPhases > lambda14Upw(0.0);
    lambda14Upw[wPhaseIdx] = (potentialDiffW14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda14Upw[nPhaseIdx] = (potentialDiffNw14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

    //compute mobilities of face 2
    Dune::FieldVector < Scalar, numPhases > lambda32Upw(0.0);
    lambda32Upw[wPhaseIdx] = (potentialDiffW32 >= 0) ? lambda3[wPhaseIdx] : lambda2[wPhaseIdx];
    lambda32Upw[nPhaseIdx] = (potentialDiffNw32 >= 0) ? lambda3[nPhaseIdx] : lambda2[nPhaseIdx];

    //compute mobilities of face 3
    Dune::FieldVector < Scalar, numPhases > lambda34Upw(0.0);
    lambda34Upw[wPhaseIdx] = (potentialDiffW34 >= 0) ? lambda3[wPhaseIdx] : lambda4[wPhaseIdx];
    lambda34Upw[nPhaseIdx] = (potentialDiffNw34 >= 0) ? lambda3[nPhaseIdx] : lambda4[nPhaseIdx];

    for (int i = 0; i < numPhases; i++)
    {
        // evaluate parts of velocity
        DimVector vel12 = interactionVolume.getNormal(0, 0);
        DimVector vel14 = interactionVolume.getNormal(3, 0);
        DimVector vel23 = interactionVolume.getNormal(1, 0);
        DimVector vel21 = interactionVolume.getNormal(0, 0);
        DimVector vel34 = interactionVolume.getNormal(2, 0);
        DimVector vel32 = interactionVolume.getNormal(1, 0);
        DimVector vel41 = interactionVolume.getNormal(3, 0);
        DimVector vel43 = interactionVolume.getNormal(2, 0);

        Dune::FieldVector < Scalar, 2 * dim > flux(0);
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

        vel12 *= flux[0] / (2 * interactionVolume.getFaceArea(0, 0)); //divide by 2 because the flux is related to the half face!
        vel14 *= flux[3] / (2 * interactionVolume.getFaceArea(0, 1));
        vel23 *= flux[1] / (2 * interactionVolume.getFaceArea(1, 0));
        vel21 *= flux[0] / (2 * interactionVolume.getFaceArea(1, 1));
        vel34 *= flux[2] / (2 * interactionVolume.getFaceArea(2, 0));
        vel32 *= flux[1] / (2 * interactionVolume.getFaceArea(2, 1));
        vel41 *= flux[3] / (2 * interactionVolume.getFaceArea(3, 0));
        vel43 *= flux[2] / (2 * interactionVolume.getFaceArea(3, 1));

        if (level1 < level2)
        {
            vel12 *= 0.5;
        }
        else if (level2 < level1)
        {
            vel21 *= 0.5;
        }
        if (level2 < level3)
        {
            vel23 *= 0.5;
        }
        else if (level3 < level2)
        {
            vel32 *= 0.5;
        }
        if (level3 < level4)
        {
            vel34 *= 0.5;
        }
        else if (level4 < level3)
        {
            vel43 *= 0.5;
        }
        if (level4 < level1)
        {
            vel41 *= 0.5;
        }
        else if (level1 < level4)
        {
            vel14 *= 0.5;
        }

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

        if (innerBoundaryVolumeFaces[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
        {
            vel12 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 1)])
        {
            vel14 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
        {
            vel23 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 1)])
        {
            vel21 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
        {
            vel34 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 1)])
        {
            vel32 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
        {
            vel41 *= 2;
        }
        if (innerBoundaryVolumeFaces[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 1)])
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
void FvMpfaL2dVelocity2p<TypeTag>::calculateBoundaryInteractionVolumeVelocity(InteractionVolume& interactionVolume,
                                                                              CellData& cellData, int elemIdx)
{
        auto element = interactionVolume.getSubVolumeElement(elemIdx);

        // get global coordinate of cell centers
        const GlobalPosition& globalPos = element.geometry().center();

        //permeability vector at boundary
        DimMatrix permeability(problem_.spatialParams().intrinsicPermeability(element));

        //get mobilities of the phases
        Dune::FieldVector < Scalar, numPhases > lambda(cellData.mobility(wPhaseIdx));
        lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

        for (int fIdx = 0; fIdx < dim; fIdx++)
        {
            int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);

            if (interactionVolume.isBoundaryFace(intVolFaceIdx))
            {
                if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(pressureEqIdx))
                {
                    int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, fIdx);

                    const auto refElement = referenceElement(element);

                    const LocalPosition& localPos = refElement.position(boundaryFaceIdx, 1);

                    const GlobalPosition& globalPosFace = element.geometry().global(localPos);

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

                    Dune::FieldVector <Scalar, numPhases> lambdaBound(fluidMatrixInteraction.krw(satWBound));
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

                    Scalar potentialDiffW = (cellData.potential(wPhaseIdx) - potentialBoundW) / dist;
                    Scalar  potentialDiffNw = (cellData.potential(nPhaseIdx) - potentialBoundNw) / dist;

                    //store potentials for further calculations (saturation, ...)
                    cellData.fluxData().addUpwindPotential(wPhaseIdx, boundaryFaceIdx, potentialDiffW);
                    cellData.fluxData().addUpwindPotential(nPhaseIdx, boundaryFaceIdx, potentialDiffNw);

                    //calculated phase velocities from advective velocities -> capillary pressure velocity already added in pressure part!
                    DimVector velocityW(0);
                    DimVector velocityNw(0);

                    // calculate capillary pressure gradient
                    DimVector pressGradient = unitDistVec;
                    pressGradient *= (cellData.potential(wPhaseIdx) - potentialBoundW) / dist;
                    permeability.mv(pressGradient, velocityW);

                    pressGradient = unitDistVec;
                    pressGradient *= (cellData.potential(nPhaseIdx) - potentialBoundNw) / dist;
                    permeability.mv(pressGradient, velocityNw);

                    velocityW *= (potentialDiffW >= 0.) ? lambda[wPhaseIdx] : lambdaBound[wPhaseIdx];
                    velocityNw *= (potentialDiffNw >= 0.) ? lambda[nPhaseIdx] : lambdaBound[nPhaseIdx];

                    //velocity is calculated from two vertices of one intersection!
                    velocityW *= 0.5;
                    velocityNw *= 0.5;

                    //store velocities
                        velocityW += cellData.fluxData().velocity(wPhaseIdx, boundaryFaceIdx);
                        velocityNw += cellData.fluxData().velocity(nPhaseIdx, boundaryFaceIdx);
                        cellData.fluxData().setVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                        cellData.fluxData().setVelocity(nPhaseIdx, boundaryFaceIdx, velocityNw);
                        cellData.fluxData().setVelocityMarker(boundaryFaceIdx);
                }
                else if (interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressureEqIdx))
                {
                    int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, fIdx);

                    const auto refElement = referenceElement(element);

                    const LocalPosition& localPos = refElement.position(boundaryFaceIdx, 1);

                    const GlobalPosition& globalPosFace = element.geometry().global(localPos);

                    DimVector distVec(globalPosFace - globalPos);
                    Scalar dist = distVec.two_norm();
                    DimVector unitDistVec(distVec);
                    unitDistVec /= dist;

                    // get neumann boundary value
                    PrimaryVariables boundValues(interactionVolume.getNeumannValues(intVolFaceIdx));

                    boundValues[wPhaseIdx] /= density_[wPhaseIdx];
                    boundValues[nPhaseIdx] /= density_[nPhaseIdx];

                    DimVector velocityW(unitDistVec);
                    DimVector velocityNw(unitDistVec);

                    velocityW *= boundValues[wPhaseIdx] / (2 * interactionVolume.getFaceArea(elemIdx, fIdx));
                    velocityNw *= boundValues[nPhaseIdx]
                            / (2 * interactionVolume.getFaceArea(elemIdx, fIdx));

                    //store potentials for further calculations (saturation, ...)
                    cellData.fluxData().addUpwindPotential(wPhaseIdx, boundaryFaceIdx, boundValues[wPhaseIdx]);
                    cellData.fluxData().addUpwindPotential(nPhaseIdx, boundaryFaceIdx, boundValues[nPhaseIdx]);

                    //store velocities
                    velocityW += cellData.fluxData().velocity(wPhaseIdx, boundaryFaceIdx);
                    velocityNw += cellData.fluxData().velocity(nPhaseIdx, boundaryFaceIdx);
                    cellData.fluxData().setVelocity(wPhaseIdx, boundaryFaceIdx, velocityW);
                    cellData.fluxData().setVelocity(nPhaseIdx, boundaryFaceIdx, velocityNw);
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

} // end namespace Dumux
#endif
