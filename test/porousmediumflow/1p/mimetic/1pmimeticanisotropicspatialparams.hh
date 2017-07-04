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
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1PMIMETICANISOTROPIC_SPATIALPARAMS_HH
#define DUMUX_1PMIMETICANISOTROPIC_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dumux
{

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(OnePMimeticTestSpatialParams);
}


/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePMimeticTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexSet = typename GridView::IndexSet;
    using ScalarVector = std::vector<Scalar>;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

    //quadrature Rule
    typedef Dune::QuadratureRule<Scalar, dim> Quad;
    typedef typename Quad::iterator QuadIterator;
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    OnePMimeticTestSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {
        perm_[0][0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, k1);
        perm_[1][1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, k2);
        perm_[1][0] = perm_[0][1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, k12);

        testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                int,
                                                Problem,
                                                TestCase);

    }

        /*!
         * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
         *
         * \param element The element
         * \param scv The sub control volume
         * \param elemSol The element solution vector
         * \return the intrinsic permeability
         */
    PermeabilityType permeability(const Element& element,
                                      const SubControlVolume& scv,
                                      const ElementSolutionVector& elemSol) const
        {
        //const GlobalPosition globalPos = element.geometry().center();
        PermeabilityType perm(0);
        perm[0][0] = 0.0;
        perm[1][1] = 0.0;
        perm[1][0] = perm[0][1] = 0.0;

        if(testCase_ != 1)
        {

            Scalar alpha_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                    Scalar,
                                                    Problem,
                                                    Alpha);

            //get the Gaussian quadrature rule for intervals
            const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());
            const Quad &quad = Dune::QuadratureRules<Scalar, dim>::rule(referenceElement.type(), 5);
            const QuadIterator qend = quad.end();

            Scalar volume = element.geometry().volume();

            QuadIterator qp = quad.begin();
            for(; qp != qend; ++qp)
            {
                GlobalPosition globalPos = element.geometry().global(qp->position());
                Scalar x = globalPos[0];
                Scalar y = globalPos[1];

                Scalar integrationElement = element.geometry().integrationElement(qp->position());

                perm[0][0] += (y*y + alpha_*x*x)/(x*x + y*y)*qp->weight()*integrationElement;
                perm[0][1] += (alpha_ - 1.0)*x*y/(x*x + y*y)*qp->weight()*integrationElement;
                perm[1][0] += (alpha_ - 1.0)*x*y/(x*x + y*y)*qp->weight()*integrationElement;
                perm[1][1] += (x*x + alpha_*y*y)/(x*x + y*y)*qp->weight()*integrationElement;
            }

            perm[0][0] /= volume;
            perm[0][1] /= volume;
            perm[1][0] /= volume;
            perm[1][1] /= volume;
        }
        else
        {
            perm = perm_;
        }

        return perm;
    }

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    PermeabilityType perm_;
    unsigned int testCase_;
};
} // end namespace
#endif
