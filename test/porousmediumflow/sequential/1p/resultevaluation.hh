// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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

#ifndef DUMUX_BENCHMARKRESULT_HH
#define DUMUX_BENCHMARKRESULT_HH

/*!
 * \file
 *
 * \ingroup IMPETtests
 * \brief Calculate errors for the diffusion test problem.
 */

#include <dumux/porousmediumflow/sequential/onemodelproblem.hh>

namespace Dumux
{
/*!
 * \brief Calculate errors for the diffusion test problem.
 */
struct ResultEvaluation
{
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

    /*!
     * \brief Calculate errors for the diffusion test problem.
     *
     * \param gridView the GridView for which the result should be evaluated
     * \param problem the Problem at hand
     * \param consecutiveNumbering indicates the order in which the
     * velocities are stored in the flux data
     */
    template<class GridView, class Problem>
    void evaluate(const GridView& gridView,
            Problem& problem, bool consecutiveNumbering = false)
    {
        using Grid = typename GridView::Grid;
        using Scalar = typename Grid::ctype;
        enum {dim=Grid::dimension};
        using Element = typename Grid::template Codim<0>::Entity;
        using Geometry = typename Element::Geometry;
        using JacobianInverseTransposed = typename Geometry::JacobianInverseTransposed;

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
        for (const auto& element : elements(gridView))
        {
            // element geometry
            const Geometry& geometry = element.geometry();
            const Dune::FieldVector<Scalar,dim>& local = referenceElement(geometry).position(0, 0);
            Dune::FieldVector<Scalar,dim> globalPos = geometry.global(local);

            Scalar volume = geometry.volume();

            int elemIdx = problem.variables().index(element);

            Scalar approxPressure = problem.variables().cellData(elemIdx).globalPressure();
            Scalar exactPressure = problem.exact(globalPos);

            numerator += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);
            denominator += volume*exactPressure*exactPressure;

            // update uMin and uMax
            using std::min;
            using std::max;
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
                double approximateFlux = problem.variables().cellData(elemIdx).fluxData().velocityTotal(isIdx)*unitOuterNormal;

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
                fluxVector[fIdx] = approximateFlux;

                if (!intersection.neighbor()) {
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

            // calculate velocity on reference element as the Raviart-Thomas-0
            // interpolant of the fluxes
            Dune::FieldVector<Scalar, dim> refVelocity;
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
                for (int i = 0; i < dim; i++)
                    refVelocity[i] = 0.5 * (fluxVector[2*i + 1] - fluxVector[2*i]);
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
        }

        using std::sqrt;
        using std::abs;
        relativeL2Error = sqrt(numerator/denominator);
        ergrad = sqrt(numeratorGrad/denominatorGrad);
        ervell2 = sqrt(numeratorFlux/denominatorFlux);
        sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf;
        errflx0 = abs((flux0 + exactflux0)/exactflux0);
        errflx1 = abs((flux1 + exactflux1)/exactflux1);
        errfly0 = abs((fluy0 + exactfluy0)/exactfluy0);
        errfly1 = abs((fluy1 + exactfluy1)/exactfluy1);
    }
};



}

#endif
