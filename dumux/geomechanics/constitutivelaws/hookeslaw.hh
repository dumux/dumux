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
 * \brief This file contains the data which is required to calculate
 *        the mechanic stresses according to Hooke's law.
 */
#ifndef DUMUX_GEOMECHANICS_HOOKES_LAW_HH
#define DUMUX_GEOMECHANICS_HOOKES_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>


namespace Dumux
{

/*!
 * \ingroup FemHookesLaw
 * \brief Evaluates the stress tensor according to hookes law for finite element models.
 */
template <class TypeTag>
class HookesLaw
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static_assert(dim == dimWorld, "Hooke's Law only works for dim = dimWorld");
    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    using LameParams = Dune::FieldVector<Scalar, 2>;

public:

    static StressTensor stressTensor(const Element& element,
                                     const IpData& ipData,
                                     const SecondaryVariables& secVars,
                                     const ElementSolution& elemSol,
                                     const LameParams& lameParams)
    {
        // evaluate displacement gradient
        StressTensor gradU(0.0);
        for (int dir = 0; dir < dim; ++dir)
            for (unsigned int i = 0; i < elemSol.size(); ++i){
                gradU[dir].axpy(elemSol[i][Indices::u(dir)], ipData.shapeGradients(i));
           //     std::cout << "i: " << i << std::endl;
           //     std::cout << "Indices::u(dir): " << Indices::u(dir) << std::endl;
           //     std::cout << "ipData.shapeValues(i): " << ipData.shapeValues(i) << std::endl;
           //     std::cout << "ipData.shapeGradients(i): " << ipData.shapeGradients(i) << std::endl;
           //     std::cout << "elemSol[i][Indices::u(dir)]: " << elemSol[i][Indices::u(dir)] << std::endl;
           //     printmatrix(std::cout, gradU, "gradU: ", "");
           //     printvector(std::cout, elemSol, "elemSol: ", "");

            }

        // evaluate strain tensor
        StressTensor epsilon;
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dimWorld; ++j)
                epsilon[i][j] = 0.5*(gradU[i][j] + gradU[j][i]);

        // trace of strain tensor
        Scalar trace(0.0);
        for (int i = 0; i < dim; ++i)
            trace += epsilon[i][i];

        // calculate sigma
        StressTensor sigma(0.0);
        for (int i = 0; i < dim; ++i)
        {
            sigma[i][i] = lameParams[0]*trace;
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] += 2.0*lameParams[1]*epsilon[i][j];
        }

        return sigma;
    }
};

} // end namespace

#endif
