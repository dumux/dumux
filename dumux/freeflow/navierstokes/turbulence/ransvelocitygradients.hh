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
 * \ingroup RANSModel
 * \copydoc Dumux::RANSVelocityGradients
 */
#ifndef DUMUX_RANS_VELOCITYGRADIENTS_HH
#define DUMUX_RANS_VELOCITYGRADIENTS_HH

#include <dune/common/fmatrix.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/problem.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Helper class for calculating the velocity gradient dependent terms for the Rans models
 */
class RansVelocityGradients
{
public:

    template <class Problem, class Element, class SubControlVolume>
    static auto velocityGradientMatrix(const Problem& problem,
                                        const Element& element,
                                        const SubControlVolume& scv)
    { return problem.ccVelocityGradients(element, scv); }

    template <class Problem, class Element, class SubControlVolume>
    static auto stressTensorScalarProduct(const Problem& problem,
                                          const Element& element,
                                          const SubControlVolume& scv)
    {
        using Scalar = typename Element::Geometry::GlobalCoordinate::value_type;
        static constexpr auto dim = Element::Geometry::GlobalCoordinate::dimension;
        using DimWorldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

        Scalar stressTensorScalarProduct = 0.0;
        const DimWorldMatrix stressTensor = stressTensor_(problem, element, scv);
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
            {
                stressTensorScalarProduct += stressTensor[dimIdx][velIdx] * stressTensor[dimIdx][velIdx];
            }
        }
        return stressTensorScalarProduct;
    }

    template <class Problem, class Element, class SubControlVolume>
    static auto vorticityTensorScalarProduct(const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        using Scalar = typename Element::Geometry::GlobalCoordinate::value_type;
        static constexpr auto dim = Element::Geometry::GlobalCoordinate::dimension;
        using DimWorldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

        Scalar vorticityTensorScalarProduct = 0.0;
        const DimWorldMatrix vorticityTensor = vorticityTensor_(problem, element, scv);
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
            {
                vorticityTensorScalarProduct += vorticityTensor[dimIdx][velIdx] * vorticityTensor[dimIdx][velIdx];
            }
        }
        return vorticityTensorScalarProduct;
    }

private:

    template <class Problem, class Element, class SubControlVolume>
    static auto stressTensor_(const Problem& problem,
                              const Element& element,
                              const SubControlVolume& scv)
    {
        using Scalar = typename Element::Geometry::GlobalCoordinate::value_type;
        static constexpr auto dim = Element::Geometry::GlobalCoordinate::dimension;
        using DimWorldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

        DimWorldMatrix stressTensor(0.0);
        const DimWorldMatrix velocityGradients = velocityGradientMatrix(problem, element, scv);
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
            {
                stressTensor[dimIdx][velIdx] = 0.5 * velocityGradients[dimIdx][velIdx]
                                             + 0.5 * velocityGradients[velIdx][dimIdx];
            }
        }
        return stressTensor;
    }

    template <class Problem, class Element, class SubControlVolume>
    static auto vorticityTensor_(const Problem& problem,
                                 const Element& element,
                                 const SubControlVolume& scv)
    {
        using Scalar = typename Element::Geometry::GlobalCoordinate::value_type;
        static constexpr auto dim = Element::Geometry::GlobalCoordinate::dimension;
        using DimWorldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

        DimWorldMatrix vorticityTensor(0.0);
        const DimWorldMatrix velocityGradients = velocityGradientMatrix(problem, element, scv);
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
            {
                vorticityTensor[dimIdx][velIdx] = 0.5 * velocityGradients[dimIdx][velIdx]
                                                - 0.5 * velocityGradients[velIdx][dimIdx];
            }
        }
        return vorticityTensor;
    }

};

} // end namespace Dumux

#endif
