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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Base class for sequential 2p2c compositional problems.
 */
#ifndef DUMUX_IMPETPROBLEM_2P2C_HH
#define DUMUX_IMPETPROBLEM_2P2C_HH

#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>
#include <dumux/porousmediumflow/sequential/variableclass.hh>
#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>


namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief  Base class for all compositional 2-phase problems which use an impet algorithm
 *
 * Extends IMPESProblem2P by the compositional the boundary formulation and initial conditions.
 * These can be specified via a feed mass fractions \f$ Z^k \f$ or a saturation, specified by
 * the appropriate flag.
 */
template<class TypeTag>
class IMPETProblem2P2C : public IMPESProblem2P<TypeTag>
{
    using ParentType = IMPESProblem2P<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;

    // material properties
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    enum
    {
        adaptiveGrid = getPropValue<TypeTag, Properties::AdaptiveGrid>()
    };
    enum {
        dimWorld = Grid::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The standard constructor
     *
     * \param timeManager The time manager
     * \param grid The grid
     */
    IMPETProblem2P2C(TimeManager& timeManager, Grid& grid)
    : IMPETProblem2P2C(timeManager, grid, grid.leafGridView())
    {}

    IMPETProblem2P2C(TimeManager &timeManager, Grid& grid, const GridView& gridView)
        : ParentType(timeManager, grid, gridView)
    { }
    /*!
     * \brief The constructor for given spatialParams
     *
     * This constructor uses a predefined SpatialParams object that was created (e.g. in
     * the problem) and does not create one in the base class.
     *
     * \param timeManager The time manager
     * \param grid The grid
     * \param spatialParams SpatialParams instantiation
     */
    IMPETProblem2P2C(TimeManager &timeManager, Grid& grid, SpatialParams &spatialParams)
        : ParentType(timeManager, grid, spatialParams)
    { }

    /*!
     * \brief Called by TimeManager just before the time
     *        integration.
     *
     *        In compositional/compressible models, the secondary variables
     *        should be updated after each time step. Hence, another update
     *        is only necessary if the grid is adapted.
     */
    void preTimeStep()
    {
        // if adaptivity is used, this method adapts the grid.
        // if it is not used, this method does nothing.
        if (this->adaptiveGrid)
        {
            this->gridAdapt().adaptGrid();

            if(this->gridAdapt().wasAdapted())
                asImp_().pressureModel().updateMaterialLaws(false);
        }
    }
    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     *
     *        In compositional/compressible models, the secondary variables
     *        should be updated for output and the next time step via
     *        updateMaterialLaws. The boolean indicates that it is
     *        a call from postTimeStep().
     */
    void postTimeStep()
    {
        ParentType::postTimeStep();
        asImp_().pressureModel().updateMaterialLaws(true);
    }
    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Saturation initial condition (dimensionless)
     *
     * The problem is initialized with the following saturation. Both
     * phases are assumed to contain an equilibrium concentration of the
     * correspondingly other component.
     * \param element The element.
     */
    Scalar initSat(const Element& element) const
    {
        return asImp_().initSatAtPos(element.geometry().center());
    }

    /*!
     * \brief Saturation initial condition (dimensionless) at given position
     *
     * Has to be provided if initSat() is not used in the specific problem.
     *  \param globalPos The global position.
     */
    Scalar initSatAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::NotImplemented, "please specify initial saturation in the problem"
                                            " using an initSatAtPos() method!");
        return NAN;
    }

    /*!
     * \brief Concentration initial condition (dimensionless)
     *
     * The problem is initialized with a  feed mass fraction:
     * Mass of component 1 per total mass \f$\mathrm{[-]}\f$. This directly
     * enters the flash calucalation.
     * \param element The element.
     */
    Scalar initConcentration(const Element& element) const
    {
        return asImp_().initConcentrationAtPos(element.geometry().center());
    }

    /*!
     * \brief Concentration initial condition (dimensionless)
     *
     * Has to be provided if initConcentration() is not used in the specific problem.
     *  \param globalPos The global position.
     */
    Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::NotImplemented, "please specify initial Concentration in the problem"
                                            " using an initConcentrationAtPos() method!");
    }
    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
    //! Sets entries of the primary variable vector to zero
    void setZero(GetPropType<TypeTag, Properties::PrimaryVariables> &values, const int equation = -1) const
    {
        if (equation == Indices::pressureEqIdx)
        {
            values[Indices::pressureEqIdx] = 0.;
        }
        else if(equation == Indices::contiNEqIdx)
            values[Indices::contiNEqIdx] =0.;
        else if(equation == Indices::contiWEqIdx)
            values[Indices::contiWEqIdx] =0.;
        else if (equation == -1)
        {
            // set everything to zero
            for (unsigned int i = 0; i < values.size(); i++)
                values[i] = 0.;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "vector of primary variables can not be set properly");
    }
};

} // end namespace Dumux

#endif
