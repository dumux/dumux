// $Id$
/****************************************************************************
*   Copyright (C) 2007-2010 by Bernd Flemisch                               *
*   Institute of Hydraulic Engineering                                      *
*   University of Stuttgart, Germany                                        *
*   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
/*!
 * \file
 *
 * \ingroup IMPETtests
 * \brief Calculate errors for a FVCA5 benchmark problem.
 */
#ifndef DUMUX_BENCHMARKRESULT_HH
#define DUMUX_BENCHMARKRESULT_HH

namespace Dumux
{

/*!
 * \brief calculate errors for a FVCA5 benchmark problem
 */
struct BenchmarkResult
{
private:
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

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
    double uMin;
    double uMax;
    double flux0;
    double flux1;
    double fluy0;
    double fluy1;
    double sumf;
    double sumflux;
    double exactflux0;
    double exactflux1;
    double exactfluy0;
    double exactfluy1;
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
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename Entity::Geometry Geometry;
        typedef typename Grid::LevelGridView GV;
        typedef typename GV::IndexSet IS;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename GV::IntersectionIterator IntersectionIterator;
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> EM;
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,FaceLayout> FM;
        typedef typename Grid::ctype ct;

        enum{dim = Grid::dimension};

        const GV& gridview(grid.levelView(grid.maxLevel()));
        const IS& indexset(gridview.indexSet());
        EM elementmapper(gridview);
        FM facemapper(gridview);

        uMean = 0;
        double domainVolume = 0;
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
        {
            // get entity
            const Entity& element = *it;

            // get volume
            double volume = element.geometry().volume();

            // cell index
            int indexi = elementmapper.map(element);

            // get approximate solution value
            uMean += volume*(problem.variables().pressure())[indexi];

            // add to domainVolume
            domainVolume += volume;
        }
        uMean /= domainVolume;

        if (pureNeumann) {
            for (int i = 0; i < (int) problem.variables().pressure().size(); i++)
                (problem.variables().pressure())[i] -= uMean;
        }


        uMin = 1e100;
        uMax = -1e100;
        flux0 = 0;
        flux1 = 0;
        fluy0 = 0;
        fluy1 = 0;
        sumf = 0;
        sumflux = 0;
        exactflux0 = 0;
        exactflux1 = 0;
        exactfluy0 = 0;
        exactfluy1 = 0;
        erflm = 0;
        ener1 = 0;
        ener2 = 0;
        double numerator = 0;
        double denominator = 0;
        double numeratorGrad = 0;
        double denominatorGrad = 0;
        double numeratorFlux = 0;
        double denominatorFlux = 0;
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
        {
            // get entity
            const Entity& element = *it;

            // element geometry
            const Geometry& geometry = element.geometry();

            // cell geometry type
            Dune::GeometryType gt = geometry.type();

            // cell center in reference element
            const Dune::FieldVector<ct,dim>&
                local = Dune::GenericReferenceElements<ct,dim>::general(gt).position(0,0);

            // get global coordinate of cell center
            Dune::FieldVector<ct,dim> global = geometry.global(local);

            // get exact solution value
            double exactValue = problem.exact(global);

            // cell index
            int indexi = elementmapper.map(element);

            // get approximate solution value
            double approximateValue = (problem.variables().pressure())[indexi];

            // update uMin and uMax
            uMin = std::min(uMin, approximateValue);
            uMax = std::max(uMax, approximateValue);

            // cell volume, assume linear map here
            double volume = geometry.volume();

            // update sumf
            sumf += volume*(problem.source(global, element, local)[0]);

            // get the absolute permeability
            Dune::FieldMatrix<double,dim,dim> K = problem.spatialParameters().K(global, element, local);

            numerator += volume*(exactValue - approximateValue)*(exactValue - approximateValue);
            denominator += volume*exactValue*exactValue;

            int i = -1;
            Dune::FieldVector<ct,2*dim> fluxVector;
            Dune::FieldVector<ct,dim> exactGradient;
            IntersectionIterator endis = gridview.iend(element);
            for (IntersectionIterator is = gridview.ibegin(element); is!=endis; ++is)
            {
                // local number of facet
                i = is->indexInInside();

                // global number of face
                int faceIndex = facemapper.template map<1>(element, i);

                Dune::FieldVector<double,dim> faceGlobal = is->geometry().center();
                double faceVol = is->geometry().volume();

                // get normal vector
                Dune::FieldVector<double,dim> unitOuterNormal = is->centerUnitOuterNormal();

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
                erflm = std::max(erflm, fabs(fluxDiff));

                numeratorFlux += volume*fluxDiff*fluxDiff;
                denominatorFlux += volume*exactFlux*exactFlux;

                // calculate the fluxes through the element faces
                exactFlux *= faceVol;
                approximateFlux *= faceVol;
                fluxVector[i] = approximateFlux;

                //if (is.boundary()) {
                if (!is->neighbor()) {
                    if (fabs(faceGlobal[1]) < 1e-6) {
                        fluy0 += approximateFlux;
                        exactfluy0 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (fabs(faceGlobal[1] - 1.0) < 1e-6) {
                        fluy1 += approximateFlux;
                        exactfluy1 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (faceGlobal[0] < 1e-6) {
                        flux0 += approximateFlux;
                        exactflux0 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                    else if (fabs(faceGlobal[0] - 1.0) < 1e-6) {
                        flux1 += approximateFlux;
                        exactflux1 += exactFlux;
                        ener2 += -approximateFlux*approximateFace;
                    }
                }
            }

            // calculate velocity on reference element
            Dune::FieldVector<ct,dim> refVelocity;
            if (geometry.corners() == 3) {
                refVelocity[0] = 1.0/3.0*(fluxVector[0] + fluxVector[2] - 2.0*fluxVector[1]);
                refVelocity[1] = 1.0/3.0*(fluxVector[0] + fluxVector[1] - 2.0*fluxVector[2]);
            }
            else {
                refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
                refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);
            }

            // get the transposed Jacobian of the element mapping
            const Dune::FieldMatrix<ct,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(local);
            Dune::FieldMatrix<ct,dim,dim> jacobianT(jacobianInv);
            jacobianT.invert();

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector<ct,dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= geometry.integrationElement(local);

            // get the approximate gradient
            Dune::FieldVector<ct,dim> approximateGradient;
            K.solve(approximateGradient, elementVelocity);

            // get the exact gradient
            exactGradient = problem.exactGrad(global);

            // the difference between exact and approximate gradient
            Dune::FieldVector<ct,dim> gradDiff(exactGradient);
            gradDiff += approximateGradient;

            // add to energy
            ener1 += volume*(approximateGradient*elementVelocity);

            numeratorGrad += volume*(gradDiff*gradDiff);
            denominatorGrad += volume*(exactGradient*exactGradient);
        }

        relativeL2Error = sqrt(numerator/denominator);
        ergrad = sqrt(numeratorGrad/denominatorGrad);
        ervell2 = sqrt(numeratorFlux/denominatorFlux);
        sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf;
        errflx0 = fabs((flux0 + exactflux0)/exactflux0);
        errflx1 = fabs((flux1 + exactflux1)/exactflux1);
        errfly0 = fabs((fluy0 + exactfluy0)/exactfluy0);
        errfly1 = fabs((fluy1 + exactfluy1)/exactfluy1);
        eren = fabs(ener1 - ener2)/std::max(ener1, ener2);

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
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

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
    double uMin;
    double uMax;
    double flux0;
    double flux1;
    double fluy0;
    double fluy1;
    double sumf;
    double sumflux;
    double exactflux0;
    double exactflux1;
    double exactfluy0;
    double exactfluy1;
    double errflx0;
    double errflx1;
    double errfly0;
    double errfly1;
    double erflm;
    double ener1;

    template<class GridView, class Problem, class SolVector, class VelVector>
    void evaluate(const GridView& gridView,
            Problem& problem,
            const SolVector& solution,
            const VelVector& velocity, bool consecutiveNumbering = false)
    {
        typedef typename GridView::Grid Grid;
        typedef typename Grid::ctype Scalar;
        enum {dim=Grid::dimension};
        typedef typename Grid::template Codim<0>::Entity Element;
        typedef typename Element::Geometry Geometry;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
        typedef typename GridView::IntersectionIterator IntersectionIterator;
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;

        ElementMapper elementMapper(gridView);

        uMin = 1e100;
        uMax = -1e100;
        flux0 = 0;
        flux1 = 0;
        fluy0 = 0;
        fluy1 = 0;
        sumf = 0;
        sumflux = 0;
        exactflux0 = 0;
        exactflux1 = 0;
        exactfluy0 = 0;
        exactfluy1 = 0;
        erflm = 0;
        ener1 = 0;
        Scalar numerator = 0;
        Scalar denominator = 0;
        double numeratorGrad = 0;
        double denominatorGrad = 0;
        double numeratorFlux = 0;
        double denominatorFlux = 0;
        ElementIterator endEIt = gridView.template end<0>();
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != endEIt; ++eIt)
        {
            const Element& element = *eIt;

            // element geometry
            const Geometry& geometry = element.geometry();

            Dune::GeometryType geomType = geometry.type();

            const Dune::FieldVector<Scalar,dim>& local = Dune::GenericReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
            Dune::FieldVector<Scalar,dim> global = geometry.global(local);

            Scalar volume = geometry.volume();

            //int eIdx = elementMapper.map(element);
            int eIdx = problem.variables().index(element);

            Scalar approxPressure = solution[eIdx];
            Scalar exactPressure = problem.exact(global);

            numerator += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
            denominator += volume*exactPressure*exactPressure;

            // update uMin and uMax
            uMin = std::min(uMin, approxPressure);
            uMax = std::max(uMax, approxPressure);

            sumf += volume*(problem.source(global, element)[0]);

            // get the absolute permeability
            Dune::FieldMatrix<double,dim,dim> K = problem.spatialParameters().intrinsicPermeability(global, element);

            int isIdx = -1;
            Dune::FieldVector<Scalar,2*dim> fluxVector;
            Dune::FieldVector<Scalar,dim> exactGradient;
            IntersectionIterator endis = gridView.iend(element);
            for (IntersectionIterator is = gridView.ibegin(element); is!=endis; ++is)
            {
                // local number of facet
                int faceIdx = is->indexInInside();

                if (consecutiveNumbering)
                    isIdx++;
                else
                    isIdx = faceIdx;

                Dune::FieldVector<double,dim> faceGlobal = is->geometry().center();
                double faceVol = is->geometry().volume();

                // get normal vector
                Dune::FieldVector<double,dim> unitOuterNormal = is->centerUnitOuterNormal();

                // get the exact gradient
                exactGradient = problem.exactGrad(faceGlobal);

                // get the negative exact velocity
                Dune::FieldVector<double,dim> KGrad(0);
                K.umv(exactGradient, KGrad);

                // calculate the exact normal velocity
                double exactFlux = KGrad*unitOuterNormal;

                // get the approximate normalvelocity
                double approximateFlux = (velocity[eIdx][isIdx]*unitOuterNormal);

                // calculate the difference in the normal velocity
                double fluxDiff = exactFlux + approximateFlux;

                // update mean value error
                erflm = std::max(erflm, fabs(fluxDiff));

                numeratorFlux += volume*fluxDiff*fluxDiff;
                denominatorFlux += volume*exactFlux*exactFlux;

                // calculate the fluxes through the element faces
                exactFlux *= faceVol;
                approximateFlux *= faceVol;
                fluxVector[faceIdx] = approximateFlux;

                if (!is->neighbor()) {
                    if (fabs(faceGlobal[1]) < 1e-6) {
                        fluy0 += approximateFlux;
                        exactfluy0 += exactFlux;
                    }
                    else if (fabs(faceGlobal[1] - 1.0) < 1e-6) {
                        fluy1 += approximateFlux;
                        exactfluy1 += exactFlux;
                    }
                    else if (faceGlobal[0] < 1e-6) {
                        flux0 += approximateFlux;
                        exactflux0 += exactFlux;
                    }
                    else if (fabs(faceGlobal[0] - 1.0) < 1e-6) {
                        flux1 += approximateFlux;
                        exactflux1 += exactFlux;
                    }
                }
            }

            // calculate velocity on reference element
            Dune::FieldVector<Scalar,dim> refVelocity;
            if (geometry.corners() == 3) {
                refVelocity[0] = 1.0/3.0*(fluxVector[0] + fluxVector[2] - 2.0*fluxVector[1]);
                refVelocity[1] = 1.0/3.0*(fluxVector[0] + fluxVector[1] - 2.0*fluxVector[2]);
            }
            else {
                refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
                refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);
            }

            // get the transposed Jacobian of the element mapping
            const Dune::FieldMatrix<Scalar,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(local);
            Dune::FieldMatrix<Scalar,dim,dim> jacobianT(jacobianInv);
            jacobianT.invert();
            //Dune::FieldMatrix<Scalar,dim,dim> jacobianT(0);

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector<Scalar,dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= geometry.integrationElement(local);

            // get the approximate gradient
            Dune::FieldVector<Scalar,dim> approximateGradient;
            K.solve(approximateGradient, elementVelocity);

            // get the exact gradient
            exactGradient = problem.exactGrad(global);

            // the difference between exact and approximate gradient
            Dune::FieldVector<Scalar,dim> gradDiff(exactGradient);
            gradDiff += approximateGradient;

            // add to energy
            ener1 += volume*(approximateGradient*elementVelocity);

            numeratorGrad += volume*(gradDiff*gradDiff);
            denominatorGrad += volume*(exactGradient*exactGradient);
        }

        relativeL2Error = sqrt(numerator/denominator);
        ergrad = sqrt(numeratorGrad/denominatorGrad);
        ervell2 = sqrt(numeratorFlux/denominatorFlux);
        sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf;
        errflx0 = fabs((flux0 + exactflux0)/exactflux0);
        errflx1 = fabs((flux1 + exactflux1)/exactflux1);
        errfly0 = fabs((fluy0 + exactfluy0)/exactfluy0);
        errfly1 = fabs((fluy1 + exactfluy1)/exactfluy1);

        return;
    }


    template<class GridView, class Problem, class SolVector, class VelVector>
    void evaluateCP(const GridView& gridView, Problem& problem,
            const SolVector& solution, const VelVector& velocity, bool switchNormals = false)
    {
        typedef typename GridView::Grid Grid;
        typedef typename Grid::ctype Scalar;
        enum {dim=Grid::dimension, maxIntersections = 12};
        typedef typename Grid::template Codim<0>::Entity Element;
        typedef typename Element::Geometry Geometry;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
        typedef typename GridView::IntersectionIterator IntersectionIterator;
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;

        ElementMapper elementMapper(gridView);

        uMin = 1e100;
        uMax = -1e100;
        sumf = 0;
        sumflux = 0;
        erflm = 0;
        Scalar maxFlux = -1;
        Scalar numerator = 0;
        Scalar denominator = 0;
        ElementIterator endEIt = gridView.template end<0>();
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != endEIt; ++eIt)
        {
            const Element& element = *eIt;

            // element geometry
            const Geometry& geometry = element.geometry();

            Dune::GeometryType geomType = geometry.type();

            const Dune::FieldVector<Scalar,dim>& local = Dune::GenericReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
            Dune::FieldVector<Scalar,dim> global = geometry.global(local);

            Scalar volume = geometry.integrationElement(local)
                    *Dune::GenericReferenceElements<Scalar,dim>::general(geomType).volume();

            int eIdx = elementMapper.map(element);

            Scalar approxPressure = solution[eIdx];
            Scalar exactPressure = problem.exact(global);
            //std::cout << global << ": p =" << exactPressure << ", p_h = " << approxPressure << std::endl;
            numerator += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
            denominator += volume*exactPressure*exactPressure;

            // update uMin and uMax
            uMin = std::min(uMin, approxPressure);
            uMax = std::max(uMax, approxPressure);

            sumf += volume*problem.source(global, element, local)[0];

            // get the absolute permeability
            Dune::FieldMatrix<Scalar,dim,dim> K = problem.spatialParameters().K(global, element, local);

            int i = -1;
            Dune::FieldVector<Scalar,dim> exactGradient;
            IntersectionIterator endis = gridView.iend(element);
            for (IntersectionIterator is = gridView.ibegin(element); is!=endis; ++is)
            {
                // get geometry type of face
                Dune::GeometryType gtf = is->geometryInInside().type();

                // local number of facet
                i++;

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>& faceLocalNm1 = Dune::GenericReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);

                // center of face in global coordinates
                Dune::FieldVector<Scalar,dim> faceGlobal = is->geometry().global(faceLocalNm1);

                // get normal vector
                Dune::FieldVector<Scalar,dim> unitOuterNormal = is->unitOuterNormal(faceLocalNm1);

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
                Scalar approximateFlux = velocity[eIdx][i]*unitOuterNormal;

                // calculate the difference in the normal velocity
                Scalar fluxDiff = exactFlux + approximateFlux;

                //                if (std::abs(fluxDiff) > 1e-16)
                //                {
                //                    std::cout << "faceGlobal " << faceGlobal << ": exact " << exactFlux << ", approximate " << approximateFlux << std::endl;
                //                }

                // update mean value error
                erflm = std::max(erflm, fabs(fluxDiff));

                maxFlux = std::max(maxFlux, std::abs(exactFlux));

                sumflux += approximateFlux;
            }
        }

        relativeL2Error = sqrt(numerator/denominator);
        sumflux -= sumf;
        erflm /= std::max(1e-6, maxFlux);

        return;
    }
};



}

#endif
