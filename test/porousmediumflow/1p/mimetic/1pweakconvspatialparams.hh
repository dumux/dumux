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
#ifndef DUMUX_ONEPWEAKCONV_SPATIALPARAMS_HH
#define DUMUX_ONEPWEAKCONV_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

template<class TypeTag>
class OnePWeakConvSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
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

    OnePWeakConvSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {
        unsigned int testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                int,
                                                Problem,
                                                TestCase);
        Scalar k1 = 1.0;
        Scalar k2 = 10.0;
        Scalar k4 = 100.0;
        Scalar k3 = 0.0;

        if(testCase_ == 1)
        {
            k3 = 10.0;
        }
        else if(testCase_ == 2)
        {
            k3 = 1000.0;
        }
        else if(testCase_ == 3)
        {
            k3 = 100000.0;
        }

        bool useScalar = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams, useScalar);

        DimWorldMatrix I(0);
        I[0][0] = 1;
        I[1][1] = 1;

        perm1_ = perm2_ = perm3_ = perm4_ = I;

        if(useScalar)
        {
            perm1_*= k1;
            perm2_*= k2;
            perm3_*= k3;
            perm4_*= k4;
        }else{
            perm1_*= k1;
            perm4_*= k4;

            perm2_[0][0] = k2;
            perm2_[1][1] = k3;
            perm2_[0][1] = perm2_[1][0] = 0.0;

            perm3_ = perm2_;
        }


        Scalar pi = 4.0*atan(1.0);
        beta_ = 7.0/16.0*pi;

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
        const GlobalPosition globalPos = element.geometry().center();
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        int region = region_(globalPos);

        if(region == 1) return perm1_;
        else if(region == 2) return perm2_;
        else if(region == 3) return perm3_;

        return perm4_;
    }

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    int region_(const GlobalPosition& globalPos) const
    {
        int region = 1;

        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        if((y - tan(beta_)*x)<=0 && y>=0) region = 1;
        else if((y - tan(beta_)*x)>0 && y>=0) region = 2;
        else if((y - tan(beta_)*x)>0 && y<0) region = 3;
        else region = 4;

        return region;
    }

private:
    PermeabilityType perm1_;
    PermeabilityType perm2_;
    PermeabilityType perm3_;
    PermeabilityType perm4_;
    Scalar beta_;
    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace

#endif
