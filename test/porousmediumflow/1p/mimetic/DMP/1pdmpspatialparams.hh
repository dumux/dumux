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
#ifndef DUMUX_1PDMP_SPATIALPARAMS_HH
#define DUMUX_1PDMP_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(OnePDMPSpatialParams);
}


/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePDMPSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
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

public:
    // export permeability type
    using PermeabilityType = DimWorldMatrix;

    OnePDMPSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {
        testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                int,
                                                Problem,
                                                TestCase);

        Scalar pi = 4.0*atan(1.0);

        thetaOne_ = 1.0*pi/6.0;

        if(testCase_ == 2)
            thetaOne_ = (2*pi/360)*67.5;

        Scalar k1_ = 1000.0;
        Scalar k2_ = 1.0;

        Scalar cost = cos(thetaOne_);
        Scalar sint = sin(thetaOne_);

        perm_[0][0] = cost*cost*k1_ + sint*sint*k2_;
        perm_[1][1] = sint*sint*k1_ + cost*cost*k2_;
        perm_[0][1] = perm_[1][0] = cost*sint*(k1_ - k2_);

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
            return perm_;
        }

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    PermeabilityType perm_;
    Scalar thetaOne_;
    unsigned int testCase_;
};
} // end namespace
#endif
