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
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of one-phase problems
 * using a fully implicit discretization method.
 */
#ifndef DUMUX_FV_SPATIAL_PARAMS_ONE_P_HH
#define DUMUX_FV_SPATIAL_PARAMS_ONE_P_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dune/common/fmatrix.hh>

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of one-phase problems
 * using a fully implicit discretization method.
 */
template<class TypeTag>
class FVSpatialParamsOneP
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Implementation = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dimWorld>;

public:
    FVSpatialParamsOneP(const Problem& problem)
    : problemPtr_(&problem)
    {}

    /*!
     * \brief Harmonic average of a discontinuous scalar field at discontinuity interface
     *        (for compatibility reasons with the function below)
     * \return the averaged scalar
     * \param T1 first scalar parameter
     * \param T2 second scalar parameter
     * \param normal The unit normal vector of the interface
     */
    Scalar harmonicMean(const Scalar T1,
                        const Scalar T2,
                        const GlobalPosition& normal) const
    { return Dumux::harmonicMean(T1, T2); }

    /*!
     * \brief Harmonic average of a discontinuous tensorial field at discontinuity interface
     * \note We do a harmonic average of the part normal to the interface (alpha*I) and
     *       an arithmetic average of the tangential part (T - alpha*I).
     * \return the averaged tensor
     * \param T1 first tensor
     * \param T2 second tensor
     * \param normal The unit normal vector of the interface
     */
    DimWorldMatrix harmonicMean(const DimWorldMatrix& T1,
                                const DimWorldMatrix& T2,
                                const GlobalPosition& normal) const
    {
        // determine nT*k*n
        GlobalPosition tmp;
        GlobalPosition tmp2;
        T1.mv(normal, tmp);
        T2.mv(normal, tmp2);
        const Scalar alpha1 = tmp*normal;
        const Scalar alpha2 = tmp2*normal;

        const Scalar alphaHarmonic = Dumux::harmonicMean(alpha1, alpha2);
        const Scalar alphaAverage = 0.5*(alpha1 + alpha2);

        DimWorldMatrix T(0.0);
        for (int i = 0; i < dimWorld; ++i)
        {
            for (int j = 0; j < dimWorld; ++j)
            {
                T[i][j] += 0.5*(T1[i][j] + T2[i][j]);
                if (i == j)
                    T[i][j] += alphaHarmonic - alphaAverage;
            }
        }

        return T;
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    decltype(auto)
    permeability(const Element& element,
                 const SubControlVolume& scv,
                 const ElementSolutionVector& elemSol) const
    {
        return asImp_().permeabilityAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \return permeability
     * \param globalPos The position of the center of the scv
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a permeability() or permeabilityAtPos() method.");
    }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        return asImp_().porosityAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \return porosity
     * \param globalPos The position of the center of the scv
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a porosityAtPos() method.");
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidHeatCapacity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolutionVector& elemSol) const
    {
        return asImp_().solidHeatCapacityAtPos(scv.center());
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a solidHeatCapacityAtPos() method.");
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        return asImp_().solidDensityAtPos(scv.center());
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a solidDensityAtPos() method.");
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    {
        return asImp_().solidThermalConductivityAtPos(scv.center());
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a solidThermalConductivityAtPos() method.");
    }

    /*!
     * \brief Function for defining the Beavers-Joseph coefficient for multidomain
     *        problems\f$\mathrm{[-]}\f$.
     *
     * \return Beavers-Joseph coefficient \f$\mathrm{[-]}\f$
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide a beaversJosephCoeffAtPos() method.");
    }

    /*!
     * \brief Apply the Forchheimer coefficient for inertial forces
     *        calculation.
     *
     *        Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90 \cite ward1964 .
     *        Actually the Forchheimer coefficient is also a function of the dimensions of the
     *        porous medium. Taking it as a constant is only a first approximation
     *        (Nield, Bejan, Convection in porous media, 2006, p. 10 \cite nield2006 )
     *
     * \param scv The sub-control volume face where the
     *           intrinsic velocity ought to be calculated.
     */
    Scalar forchCoeff(const SubControlVolumeFace &scvf) const
    {
        static Scalar forchCoeff = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "SpatialParams.ForchCoeff", 0.55);
        return forchCoeff;
    }

    //! The problem we are associated with
    const Problem& problem() const
    {
        return *problemPtr_;
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Problem *problemPtr_;
};

} // namespace Dumux

#endif
