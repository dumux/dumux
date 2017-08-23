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
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1947)
 */
#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>

// solve Navier-Stokes equations
#define ENABLE_NAVIERSTOKES 1


namespace Dumux
{
template <class TypeTag>
class KovasznayTestProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<KovasznayTestProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
NEW_TYPE_TAG(KovasznayTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

SET_PROP(KovasznayTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(KovasznayTestProblem, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(KovasznayTestProblem, Problem, Dumux::KovasznayTestProblem<TypeTag> );

SET_BOOL_PROP(KovasznayTestProblem, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(KovasznayTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(KovasznayTestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(KovasznayTestProblem, ProblemEnableGravity, true);

#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(KovasznayTestProblem, EnableInertiaTerms, true);
#else
SET_BOOL_PROP(KovasznayTestProblem, EnableInertiaTerms, false);
#endif
}

/*!
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the staggered grid (Kovasznay 1947)
 * \todo doc me!
 */
template <class TypeTag>
class KovasznayTestProblem : public NavierStokesProblem<TypeTag>
{
    typedef NavierStokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    KovasznayTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView), eps_(1e-6)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        printL2Error_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                     bool,
                                                     Problem,
                                                     PrintL2Error);

        kinematicViscosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, LiquidKinematicViscosity);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);

        using CellArray = std::array<unsigned int, dimWorld>;
        const CellArray numCells = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                      CellArray,
                                                      Grid,
                                                      Cells);
        cellSizeX_ = this->bBoxMax()[0] / numCells[0];
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    void postTimeStep() const
    {
        if(printL2Error_)
        {
            const auto l2error = calculateL2Error();
            const int numCellCenterDofs = this->model().numCellCenterDofs();
            const int numFaceDofs = this->model().numFaceDofs();
            std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                    << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                    << std::scientific
                    << "L2(p) = " << l2error.first[pressureIdx] << " / " << l2error.second[pressureIdx]
                    << ", L2(vx) = " << l2error.first[velocityXIdx] << " / " << l2error.second[velocityXIdx]
                    << ", L2(vy) = " << l2error.first[velocityYIdx] << " / " << l2error.second[velocityYIdx]
                    << std::endl;
        }
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }


    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    SourceValues sourceAtPos(const GlobalPosition &globalPos) const
    {
        return SourceValues(0.0);
    }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

        // set a fixed pressure in one cell
        if (isLowerLeftCell_(globalPos))
            values.setDirichletCell(massBalanceIdx);
        else
            values.setOutflow(massBalanceIdx);

        return values;
    }

    /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    BoundaryValues dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos);
    }

     /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    BoundaryValues analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        BoundaryValues values;
        values[pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
        values[velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
        values[velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values;
        values[pressureIdx] = 0.0;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& pressureExact = outputModule.createScalarField("pressureExact", 0);
        auto& velocityExact = outputModule.createVectorField("velocityExact", 0);

        auto& scalarFaceVelocityExact = outputModule.createFaceScalarField("scalarFaceVelocityExact");
        auto& vectorFaceVelocityExact = outputModule.createFaceVectorField("vectorFaceVelocityExact");

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = dirichletAtPos(ccDofPosition);

                GlobalPosition velocityVector(0.0);
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    auto faceDofIdx = scvf.dofIndex();
                    auto faceDofPosition = scvf.center();
                    auto dirIdx = scvf.directionIndex();
                    auto analyticalSolutionAtFace = dirichletAtPos(faceDofPosition);
                    scalarFaceVelocityExact[faceDofIdx] = analyticalSolutionAtFace[faceIdx][dirIdx];

                    GlobalPosition tmp(0.0);
                    tmp[dirIdx] = analyticalSolutionAtFace[faceIdx][dirIdx];
                    vectorFaceVelocityExact[faceDofIdx] = std::move(tmp);
                }
                pressureExact[ccDofIdx] = analyticalSolutionAtCc[pressureIdx];
                velocityExact[ccDofIdx] = analyticalSolutionAtCc[faceIdx];
            }
        }
    }

    /*!
     * \brief Calculate the L2 error between the analytical solution and the numerical approximation.
     *
     */
    auto calculateL2Error() const
    {
        BoundaryValues sumError(0.0), sumReference(0.0), l2NormAbs(0.0), l2NormRel(0.0);

        const int numFaceDofs = this->model().numFaceDofs();

        std::vector<Scalar> staggeredVolume(numFaceDofs);
        std::vector<Scalar> errorVelocity(numFaceDofs);
        std::vector<Scalar> velocityReference(numFaceDofs);
        std::vector<int> directionIndex(numFaceDofs);

        Scalar totalVolume = 0.0;

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                // treat cell-center dofs
                const auto dofIdxCellCenter = scv.dofIndex();
                const auto& posCellCenter = scv.dofPosition();
                const auto analyticalSolutionCellCenter = dirichletAtPos(posCellCenter)[cellCenterIdx];
                const auto numericalSolutionCellCenter = this->model().curSol()[cellCenterIdx][dofIdxCellCenter];
                sumError[cellCenterIdx] += squaredDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter) * scv.volume();
                sumReference[cellCenterIdx] += analyticalSolutionCellCenter * analyticalSolutionCellCenter * scv.volume();
                totalVolume += scv.volume();

                // treat face dofs
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const int dofIdxFace = scvf.dofIndex();
                    const int dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionFace = dirichletAtPos(scvf.center())[faceIdx][dirIdx];
                    const auto numericalSolutionFace = this->model().curSol()[faceIdx][dofIdxFace][momentumBalanceIdx];
                    directionIndex[dofIdxFace] = dirIdx;
                    errorVelocity[dofIdxFace] = squaredDiff_(analyticalSolutionFace, numericalSolutionFace);
                    velocityReference[dofIdxFace] = squaredDiff_(analyticalSolutionFace, 0.0);
                    const Scalar staggeredHalfVolume = 0.5 * scv.volume();
                    staggeredVolume[dofIdxFace] = staggeredVolume[dofIdxFace] + staggeredHalfVolume;
                }
            }
        }

        // get the absolute and relative discrete L2-error for cell-center dofs
        l2NormAbs[cellCenterIdx] = std::sqrt(sumError[cellCenterIdx] / totalVolume);
        l2NormRel[cellCenterIdx] = std::sqrt(sumError[cellCenterIdx] / sumReference[cellCenterIdx]);

        // get the absolute and relative discrete L2-error for face dofs
        for(int i = 0; i < numFaceDofs; ++i)
        {
            const int dirIdx = directionIndex[i];
            const auto error = errorVelocity[i];
            const auto ref = velocityReference[i];
            const auto volume = staggeredVolume[i];
            sumError[faceIdx][dirIdx] += error * volume;
            sumReference[faceIdx][dirIdx] += ref * volume;
        }

        for(int dirIdx = 0; dirIdx < dimWorld; ++dirIdx)
        {
            l2NormAbs[faceIdx][dirIdx] = std::sqrt(sumError[faceIdx][dirIdx] / totalVolume);
            l2NormRel[faceIdx][dirIdx] = std::sqrt(sumError[faceIdx][dirIdx] / sumReference[faceIdx][dirIdx]);
        }
        return std::make_pair(l2NormAbs, l2NormRel);
    }

private:
    template<class T>
    T squaredDiff_(const T& a, const T& b) const
    {
        return (a-b)*(a-b);
    }

    bool isLowerLeftCell_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < (this->bBoxMin()[0] + 0.5*cellSizeX_ + eps_);
    }

    Scalar eps_;
    Scalar cellSizeX_;
    std::string name_;

    Scalar kinematicViscosity_;
    Scalar lambda_;
    bool printL2Error_;
};
} //end namespace

#endif
