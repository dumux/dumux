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
 *
 * \ingroup IMPETtests
 * \brief Calculate errors for a FVCA6 benchmark problem.
 */
#ifndef DUMUX_BENCHMARKRESULT_HH
#define DUMUX_BENCHMARKRESULT_HH

#include <dumux/porousmediumflow/sequential/onemodelproblem.hh>

namespace Dumux
{

/*!
 * \brief calculate errors for a FVCA5 benchmark problem
 */
struct BenchmarkResult
{
private:
    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim-1;
        }
    };

public:
    double relativeL2Error;
    double ergrad;
    double ervell2;
    double uMinExact;
    double uMaxExact;
    double uMin;
    double uMax;
    double flux0;
    double flux1;
    double fluy0;
    double fluy1;
    double fluz0;
    double fluz1;
    double sumf;
    double sumflux;
    double exactflux0;
    double exactflux1;
    double exactfluy0;
    double exactfluy1;
    double exactfluz0;
    double exactfluz1;
    double errflx0;
    double errflx1;
    double errfly0;
    double errfly1;
    double erflm;
    double ener1;
    double ener2;
    double eren;
    double uMean;

    template<class Grid, class ProblemType, class SolutionType>
    void evaluate(const Grid& grid, ProblemType& problem,
                  SolutionType& solution, bool pureNeumann = false)
    {
        using Entity = typename Grid::Traits::template Codim<0>::Entity;
        using Geometry = typename Entity::Geometry;
        using GV = typename Grid::LevelGridView;
        using IS = typename GV::IndexSet;
        using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
        using ct = typename Grid::ctype;

        enum{dim = Grid::dimension};

        using JacobianInverseTransposed = typename Geometry::JacobianInverseTransposed;
        const GV& gridview(grid.levelGridView(grid.maxLevel()));
        const IS& indexset(gridview.indexSet());
        Mapper elementMapper(gridview, Dune::mcmgElementLayout());
        Mapper faceMapper(gridview, Dune::mcmgLayout(Dune::Codim<1>()));
        SolutionType& exactSol(grid.levelGridView(grid.maxLevel()));

        uMean = 0;
        double domainVolume = 0;
        for (const auto& element : elements(gridview))
        {
            // get volume
            double volume = element.geometry().volume();

            // cell index
            int indexi = elementMapper.index(element);

            // get approximate solution value
            uMean += volume*(problem.variables().cellData(indexi).globalPressure());

            // add to domainVolume
            domainVolume += volume;
        }
        uMean /= domainVolume;

        if (pureNeumann) {
            for (int i = 0; i < (int) problem.variables().pressure().size(); i++)
            {
                double press = problem.variables().cellData(i).globalPressure() - uMean;
                problem.variables().cellData(i).setGlobalPressure(press);
            }
        }


        uMinExact = 1e100;
        uMaxExact = -1e100;
        uMin = 1e100;
        uMax = -1e100;
        flux0 = 0;
        flux1 = 0;
        fluy0 = 0;
        fluy1 = 0;
        fluz0 = 0;
        fluz1 = 0;
        sumf = 0;
        sumflux = 0;
        exactflux0 = 0;
        exactflux1 = 0;
        exactfluy0 = 0;
        exactfluy1 = 0;
        exactfluz0 = 0;
        exactfluz1 = 0;
        erflm = 0;
        ener1 = 0;
        ener2 = 0;
        double numerator = 0;
        double denominator = 0;
        double numeratorGrad = 0;
        double denominatorGrad = 0;
        double numeratorFlux = 0;
        double denominatorFlux = 0;
        bool exactOutput = true;
        for (const auto& element : elements(gridview))
        {
            // element geometry
            const Geometry& geometry = element.geometry();

            // cell center in reference element
            const Dune::FieldVector<ct,dim>& local = referenceElement(geometry).position(0,0);

            // get global coordinate of cell center
            Dune::FieldVector<ct,dim> globalPos = geometry.global(local);

            // get exact solution value
            double exactValue = problem.exact(globalPos);

            // cell index
            int indexi = elementMapper.index(element);

            // evaluate exact solution vector
            exactSol[indexi] = exactValue;

            // get approximate solution value
            double approximateValue = problem.variables().cellData(indexi).globalPressure();

            // update uMinExact and uMaxExact
            using std::min;
            using std::max;
            uMinExact = min(uMinExact, exactValue);
            uMaxExact = max(uMaxExact, exactValue);

            // update uMin and uMax
            uMin = min(uMin, approximateValue);
            uMax = max(uMax, approximateValue);

            // cell volume, assume linear map here
            double volume = geometry.volume();

            // update sumf
            sumf += volume*(problem.source(globalPos, element, local)[0]);

            // get the absolute permeability
            Dune::FieldMatrix<double,dim,dim> K = problem.spatialParams().K(globalPos, element, local);

            numerator += volume*(exactValue - approximateValue)*(exactValue - approximateValue);
            denominator += volume*exactValue*exactValue;

            int i = -1;
            Dune::FieldVector<ct,2*dim> fluxVector;
            Dune::FieldVector<ct,dim> exactGradient;
            for (const auto& intersection : intersections(gridview, element))
            {
                // local number of facet
                i = intersection.indexInInside();

                // global number of face
                int faceIndex = faceMapper.template map<1>(element, i);

                Dune::FieldVector<double,dim> faceGlobal = intersection.geometry().center();
                double faceVol = intersection.geometry().volume();

                // get normal vector
                Dune::FieldVector<double,dim> unitOuterNormal = intersection.centerUnitOuterNormal();

                // get the approximate solution on the face
                double approximateFace = (*solution.pressTrace)[faceIndex];

                // get the exact gradient
                exactGradient = problem.exactGrad(faceGlobal);

                // get the negative exact velocity
                Dune::FieldVector<double,dim> KGrad(0);
                K.umv(exactGradient, KGrad);

                // calculate the exact normal velocity
                double exactFlux = KGrad*unitOuterNormal;

                // get the approximate normalvelocity
                double approximateFlux = solution.normalVelocity[indexi][i];

                // calculate the difference in the normal velocity
                double fluxDiff = exactFlux + approximateFlux;

                // update mean value error
                using std::abs;
                erflm = max(erflm, abs(fluxDiff));

                numeratorFlux += volume*fluxDiff*fluxDiff;
                denominatorFlux += volume*exactFlux*exactFlux;

                // calculate the fluxes through the element faces
                exactFlux *= faceVol;
                approximateFlux *= faceVol;
                fluxVector[i] = approximateFlux;

                //if (is.boundary()) {
                if (!intersection.neighbor()) {
                    if (abs(faceGlobal[2]) < 1e-6) {
                        fluz0 += approximateFlux;
                        exactfluz0 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (abs(faceGlobal[2] - 1.0) < 1e-6) {
                        fluz1 += approximateFlux;
                        exactfluz1 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    if (abs(faceGlobal[1]) < 1e-6) {
                        fluy0 += approximateFlux;
                        exactfluy0 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (abs(faceGlobal[1] - 1.0) < 1e-6) {
                        fluy1 += approximateFlux;
                        exactfluy1 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (faceGlobal[0] < 1e-6) {
                        flux0 += approximateFlux;
                        exactflux0 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (abs(faceGlobal[0] - 1.0) < 1e-6) {
                        flux1 += approximateFlux;
                        exactflux1 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                }
            }

            // calculate velocity on reference element as the Raviart-Thomas-0
            // interpolant of the fluxes
            Dune::FieldVector<double, dim> refVelocity;
            // simplices
            if (geometry.type().isSimplex()) {
                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                {
                    refVelocity[dimIdx] = -fluxVector[dim - 1 - dimIdx];
                    for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                    {
                        refVelocity[dimIdx] += fluxVector[fIdx]/(dim + 1);
                    }
                }
            }
            // cubes
            else if (geometry.type().isCube()){
                for (int j = 0; j < dim; j++)
                    refVelocity[j] = 0.5 * (fluxVector[2*j + 1] - fluxVector[2*j]);
            }
            // 3D prism and pyramids
            else {
                DUNE_THROW(Dune::NotImplemented, "velocity output for prism/pyramid not implemented");
            }

            // get the transposed Jacobian of the element mapping
            const JacobianInverseTransposed& jacobianInv = geometry.jacobianInverseTransposed(local);
            JacobianInverseTransposed jacobianT(jacobianInv);
            jacobianT.invert();

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector<ct,dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= geometry.integrationElement(local);

            // get the approximate gradient
            Dune::FieldVector<ct,dim> approximateGradient;
            K.solve(approximateGradient, elementVelocity);

            // get the exact gradient
            exactGradient = problem.exactGrad(globalPos);

            // the difference between exact and approximate gradient
            Dune::FieldVector<ct,dim> gradDiff(exactGradient);
            gradDiff += approximateGradient;

            // add to energy
            ener1 += volume*(approximateGradient*elementVelocity);

            numeratorGrad += volume*(gradDiff*gradDiff);
            denominatorGrad += volume*(exactGradient*exactGradient);
        }

        using std::sqrt;
        using std::abs;
        using std::max;
        relativeL2Error = sqrt(numerator/denominator);
        ergrad = sqrt(numeratorGrad/denominatorGrad);
        ervell2 = sqrt(numeratorFlux/denominatorFlux);
        sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf;
        errflx0 = abs((flux0 + exactflux0)/exactflux0);
        errflx1 = abs((flux1 + exactflux1)/exactflux1);
        errfly0 = abs((fluy0 + exactfluy0)/exactfluy0);
        errfly1 = abs((fluy1 + exactfluy1)/exactfluy1);
        eren = abs(ener1 - ener2)/max(ener1, ener2);

        // generate a VTK file of exact solution on the maximal grid level
        if (exactOutput)
        {
            Dune::VTKWriter<GV> vtkwriter(gridview);
            char fname[128];
            sprintf(fname,"%d","exactSol-numRefine", grid.maxLevel());
            vtkwriter.addCellData(exactSol,"exact pressure solution~");
            vtkwriter.write(fname,Dune::VTK::ascii);
        }

        return;
    }
};

/*!
 * \brief calculate errors for a FVCA5 benchmark problem
 */
struct ResultEvaluation
{
private:
    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim-1;
        }
    };

public:
    double relativeL2Error;
    double relativeL2ErrorIn;
    double relativeL2ErrorOut;
    double ergrad;
    double ervell2;
    double ervell2In;
    double uMinExact;
    double uMaxExact;
    double uMin;
    double uMax;
    double flux0;
    double flux1;
    double fluy0;
    double fluy1;
    double fluz0;
    double fluz1;
    double sumf;
    double sumflux;
    double exactflux0;
    double exactflux1;
    double exactfluy0;
    double exactfluy1;
    double exactfluz0;
    double exactfluz1;
    double errflx0;
    double errflx1;
    double errfly0;
    double errfly1;
    double erflm;
    double ener1;
    double hMax;

    template<class GridView, class Problem>
    void evaluate(const GridView& gridView,
            Problem& problem, bool consecutiveNumbering = false)
    {
        using Grid = typename GridView::Grid;
        using Scalar = typename Grid::ctype;
        enum {dim=Grid::dimension};
        using Element = typename Grid::template Codim<0>::Entity;
        using Geometry = typename Element::Geometry;
        using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
        using SolVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
        using JacobianInverseTransposed = typename Geometry::JacobianInverseTransposed;

        ElementMapper elementMapper(gridView, Dune::mcmgElementLayout());
        SolVector exactSol(gridView.size(0));

        uMinExact = 1e100;
        uMaxExact = -1e100;
        uMin = 1e100;
        uMax = -1e100;
        flux0 = 0;
        flux1 = 0;
        fluy0 = 0;
        fluy1 = 0;
        fluz0 = 0;
        fluz1 = 0;
        sumf = 0;
        sumflux = 0;
        exactflux0 = 0;
        exactflux1 = 0;
        exactfluy0 = 0;
        exactfluy1 = 0;
        exactfluz0 = 0;
        exactfluz1 = 0;
        erflm = 0;
        ener1 = 0;
        hMax = 0;
        Scalar numerator = 0;
        Scalar denominator = 0;
        Scalar numeratorIn = 0;
        Scalar denominatorIn = 0;
        Scalar numeratorOut = 0;
        Scalar denominatorOut = 0;
        double numeratorGrad = 0;
        double denominatorGrad = 0;
        double numeratorFlux = 0;
        double denominatorFlux = 0;
        double numeratorFluxIn = 0;
        double denominatorFluxIn = 0;
        bool exactOutput = true;
        for (const auto& element : elements(gridView))
        {
            // element geometry
            const Geometry& geometry = element.geometry();
            const Dune::FieldVector<Scalar,dim>& local = referenceElement(geometry).position(0, 0);
            Dune::FieldVector<Scalar,dim> globalPos = geometry.global(local);

            Scalar volume = geometry.volume();

            int eIdx = problem.variables().index(element);

            Scalar approxPressure = problem.variables().cellData(eIdx).globalPressure();
            Scalar exactPressure = problem.exact(globalPos);

            // evaluate exact solution vector
            exactSol[eIdx] = exactPressure;

            // output local relative error for each cell
            //std::cout << "local relative error for cell "<< eIdx << " is: "
            //          << (approxPressure - exactPressure)/exactPressure << std::endl;

            numerator += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
            denominator += volume*exactPressure*exactPressure;

            // calculate the relative error for inner cells
            bool bndCell = false;
            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                {
                    bndCell = true;
                    break;
                }
            }

            if (bndCell)
            {
                numeratorOut += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
                denominatorOut += volume*exactPressure*exactPressure;
            }
            else
            {
                numeratorIn += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
                denominatorIn += volume*exactPressure*exactPressure;
            }

            // update uMinExact and uMaxExact
            using std::min;
            using std::max;
            uMinExact = min(uMinExact, exactPressure);
            uMaxExact = max(uMaxExact, exactPressure);

            // update uMin and uMax
            uMin = min(uMin, approxPressure);
            uMax = max(uMax, approxPressure);

            typename Problem::PrimaryVariables sourceVec;
            problem.source(sourceVec, element);
            sumf += volume*sourceVec[0];

            // get the absolute permeability
            Dune::FieldMatrix<double,dim,dim> K = problem.spatialParams().intrinsicPermeability(element);

            int isIdx = -1;
            Dune::FieldVector<Scalar,2*dim> fluxVector;
            Dune::FieldVector<Scalar,dim> exactGradient;
            for (const auto& intersection : intersections(gridView, element))
            {
                // local number of facet
                int fIdx = intersection.indexInInside();

                if (consecutiveNumbering)
                    isIdx++;
                else
                    isIdx = fIdx;

                Dune::FieldVector<double,dim> faceGlobal = intersection.geometry().center();
                double faceVol = intersection.geometry().volume();

                // get normal vector
                Dune::FieldVector<double,dim> unitOuterNormal = intersection.centerUnitOuterNormal();

                // get the exact gradient
                exactGradient = problem.exactGrad(faceGlobal);

                // get the negative exact velocity
                Dune::FieldVector<double,dim> KGrad(0);
                K.umv(exactGradient, KGrad);

                // calculate the exact normal velocity
                double exactFlux = KGrad*unitOuterNormal;

                // get the approximate normalvelocity
                double approximateFlux = problem.variables().cellData(eIdx).fluxData().velocityTotal(fIdx)*unitOuterNormal;

                // calculate the difference in the normal velocity
                double fluxDiff = exactFlux + approximateFlux;

                // update mean value error
                using std::abs;
                using std::max;
                erflm = max(erflm, abs(fluxDiff));

                numeratorFlux += volume*fluxDiff*fluxDiff;
                denominatorFlux += volume*exactFlux*exactFlux;

                if (!bndCell)
                {
                    numeratorFluxIn += volume*fluxDiff*fluxDiff;
                    denominatorFluxIn += volume*exactFlux*exactFlux;
                }

                // calculate the fluxes through the element faces
                exactFlux *= faceVol;
                approximateFlux *= faceVol;
                fluxVector[fIdx] = approximateFlux;

                if (!intersection.neighbor()) {
                    if (abs(faceGlobal[2]) < 1e-6) {
                        fluz0 += approximateFlux;
                        exactfluz0 += exactFlux;
                    }
                    else if (abs(faceGlobal[2] - 1.0) < 1e-6) {
                        fluz1 += approximateFlux;
                        exactfluz1 += exactFlux;
                    }
                    if (abs(faceGlobal[1]) < 1e-6) {
                        fluy0 += approximateFlux;
                        exactfluy0 += exactFlux;
                    }
                    else if (abs(faceGlobal[1] - 1.0) < 1e-6) {
                        fluy1 += approximateFlux;
                        exactfluy1 += exactFlux;
                    }
                    else if (faceGlobal[0] < 1e-6) {
                        flux0 += approximateFlux;
                        exactflux0 += exactFlux;
                    }
                    else if (abs(faceGlobal[0] - 1.0) < 1e-6) {
                        flux1 += approximateFlux;
                        exactflux1 += exactFlux;
                    }
                }
            }

            // calculate velocity on reference element
            Dune::FieldVector<Scalar,dim> refVelocity;
            if (geometry.corners() == 8) {
                refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
                refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);
                refVelocity[2] = 0.5*(fluxVector[5] - fluxVector[4]);
            }
            else {

            }

            // get the transposed Jacobian of the element mapping
            const JacobianInverseTransposed& jacobianInv = geometry.jacobianInverseTransposed(local);
            JacobianInverseTransposed jacobianT(jacobianInv);
            jacobianT.invert();

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector<Scalar,dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= geometry.integrationElement(local);

            // get the approximate gradient
            Dune::FieldVector<Scalar,dim> approximateGradient;
            K.solve(approximateGradient, elementVelocity);

            // get the exact gradient
            exactGradient = problem.exactGrad(globalPos);

            // the difference between exact and approximate gradient
            Dune::FieldVector<Scalar,dim> gradDiff(exactGradient);
            gradDiff += approximateGradient;

            // add to energy
            ener1 += volume*(approximateGradient*elementVelocity);

            numeratorGrad += volume*(gradDiff*gradDiff);
            denominatorGrad += volume*(exactGradient*exactGradient);

            // calculate the maximum of the diagonal length of all elements on leaf grid
            for (int i = 0; i < element.geometry().corners(); ++i)
            {
                Dune::FieldVector<Scalar,dim> corner1 = element.geometry().corner(i);

                for (int j = 0; j < element.geometry().corners(); ++j)
                {
                    // get all corners of current element and compare the distances between them
                    Dune::FieldVector<Scalar,dim> corner2 = element.geometry().corner(j);

                    // distance vector between corners
                    Dune::FieldVector<Scalar,dim> distVec = corner1 - corner2;
                    Scalar dist = distVec.two_norm();

                    if (hMax < dist)
                        hMax = dist;
                }
            }
        }

        using std::sqrt;
        using std::abs;
        relativeL2Error = sqrt(numerator/denominator);
        relativeL2ErrorIn = sqrt(numeratorIn/denominatorIn);
        relativeL2ErrorOut = sqrt(numeratorOut/denominatorOut);
        ergrad = sqrt(numeratorGrad/denominatorGrad);
        ervell2 = sqrt(numeratorFlux/denominatorFlux);
        ervell2In = sqrt(numeratorFluxIn/denominatorFluxIn);
        sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf;
        errflx0 = abs((flux0 + exactflux0)/exactflux0);
        errflx1 = abs((flux1 + exactflux1)/exactflux1);
        errfly0 = abs((fluy0 + exactfluy0)/exactfluy0);
        errfly1 = abs((fluy1 + exactfluy1)/exactfluy1);

        // generate a VTK file of exact solution
        if (exactOutput)
        {
            Dune::VTKWriter<GridView> vtkwriter(gridView);
            char fname[128];
            sprintf(fname, "exactSol-numRefine%d", gridView.grid().maxLevel());
            vtkwriter.addCellData(exactSol,"exact pressure solution~");
            vtkwriter.write(fname,Dune::VTK::ascii);
        }

        return;
    }


    template<class GridView, class Problem, class SolVector, class VelVector>
    void evaluateCP(const GridView& gridView, Problem& problem,
            const SolVector& solution, const VelVector& velocity, bool switchNormals = false)
    {
        using Grid = typename GridView::Grid;
        using Scalar = typename Grid::ctype;
        enum {dim=Grid::dimension, maxIntersections = 12};
        using Element = typename Grid::template Codim<0>::Entity;
        using Geometry = typename Element::Geometry;
        using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

        ElementMapper elementMapper(gridView, Dune::mcmgElementLayout());

        uMinExact = 1e100;
        uMaxExact = -1e100;
        uMin = 1e100;
        uMax = -1e100;
        sumf = 0;
        sumflux = 0;
        erflm = 0;
        Scalar maxFlux = -1;
        Scalar numerator = 0;
        Scalar denominator = 0;
        for (const auto& element : elements(gridView))
        {
            // element geometry
            const Geometry& geometry = element.geometry();
            const Dune::FieldVector<Scalar,dim>& local = referenceElement(geometry).position(0, 0);
            Dune::FieldVector<Scalar,dim> globalPos = geometry.global(local);

            Scalar volume = geometry.integrationElement(local) * referenceElement(geometry).volume();

            int eIdx = elementMapper.index(element);

            Scalar approxPressure = solution[eIdx];
            Scalar exactPressure = problem.exact(globalPos);
            numerator += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
            denominator += volume*exactPressure*exactPressure;

            // update uMinExact and uMaxExact
            using std::min;
            using std::max;
            uMinExact = min(uMinExact, exactPressure);
            uMaxExact = max(uMaxExact, exactPressure);

            // update uMin and uMax
            uMin = min(uMin, approxPressure);
            uMax = max(uMax, approxPressure);

            sumf += volume*problem.source(globalPos, element, local)[0];

            // get the absolute permeability
            Dune::FieldMatrix<Scalar,dim,dim> K = problem.spatialParams().K(globalPos, element, local);

            int i = -1;
            Dune::FieldVector<Scalar,dim> exactGradient;
            for (const auto& intersection : intersections(gridView, element))
            {
                // local number of facet
                i++;

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>& faceLocalNm1 = referenceElement(intersection.geometryInInside()).position(0,0);

                // center of face in global coordinates
                Dune::FieldVector<Scalar,dim> faceGlobal = intersection.geometry().global(faceLocalNm1);

                // get normal vector
                Dune::FieldVector<Scalar,dim> unitOuterNormal = intersection.unitOuterNormal(faceLocalNm1);

                if (switchNormals)
                    unitOuterNormal *= -1.0;

                // get the exact gradient
                exactGradient = problem.exactGrad(faceGlobal);
                exactGradient.axpy(-problem.wettingPhase().density(), problem.gravity());

                // get the negative exact velocity
                Dune::FieldVector<Scalar,dim> KGrad(0);
                K.umv(exactGradient, KGrad);

                // calculate the exact normal velocity
                Scalar exactFlux = KGrad*unitOuterNormal;

                // get the approximate normalvelocity
                Scalar approximateFlux = problem.cellData(eIdx).fluxData().velocityTotal(i)*unitOuterNormal;

                // calculate the difference in the normal velocity
                Scalar fluxDiff = exactFlux + approximateFlux;

                using std::abs;
                //                if (abs(fluxDiff) > 1e-16)
                //                {
                //                    std::cout << "faceGlobal " << faceGlobal << ": exact "
                //                              << exactFlux << ", approximate " << approximateFlux << std::endl;
                //                }

                // update mean value error
                erflm = max(erflm, abs(fluxDiff));

                maxFlux = max(maxFlux, abs(exactFlux));

                sumflux += approximateFlux;
            }
        }

        using std::sqrt;
        using std::max;
        relativeL2Error = sqrt(numerator/denominator);
        sumflux -= sumf;
        erflm /= max(1e-6, maxFlux);

        return;
    }
};

}

#endif
