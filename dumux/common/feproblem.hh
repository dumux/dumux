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
 * \ingroup Common
 * \brief Base class for all finite element problems
 */
#ifndef DUMUX_COMMON_FE_PROBLEM_HH
#define DUMUX_COMMON_FE_PROBLEM_HH

#include <memory>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>

namespace Dumux {

//! \todo TODO: Put generic problem traits somewhere central!
template<class Problem>
struct ProblemTraits;

/*!
 * \ingroup Common
 * \brief Base class for all problems solved on the basis of a finite-element scheme.
 * \tparam GridGeometry The grid geometry class
 * \tparam PrimaryVariables The type used for primary variables
 * \tparam Impl The actual problem implementation
 *
 * \note All quantities (regarding the units) are specified assuming a
 *       three-dimensional world. Problems discretized using 2D grids
 *       are assumed to be extruded by \f$1 m\f$ and 1D grids are assumed
 *       to have a cross section of \f$1m \times 1m\f$.
 *
 * \todo TODO: Implement support for point sources.
 * \todo TODO: Implement support for internal Dirichlet constraints.
 * \todo TODO: Add interfaces for analytic source term derivatives.
 */
template<class Impl>
class FEProblem
{
    using GridGeometry = typename ProblemTraits<Impl>::GridGeometry;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Intersection = typename GridView::Intersection;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = typename ProblemTraits<Impl>::BoundaryTypes;
    using NumEqVector = typename ProblemTraits<Impl>::NumEqVector;
    using PrimaryVariables = typename ProblemTraits<Impl>::PrimaryVariables;
    using Scalar = typename ProblemTraits<Impl>::Scalar;

    static constexpr int dim = GridView::dimension;

public:
    //! export traits of this problem
    using Traits = ProblemTraits<Impl>;

    /*!
     * \brief Constructor
     * \param gridGeometry The grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first
     */
    FEProblem(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , paramGroup_(paramGroup)
    {
        // set a default name for the problem
        problemName_ = getParamFromGroup<std::string>(paramGroup, "Problem.Name");
        static_assert(!Impl::enableInternalDirichletConstraints(),
                      "Internal Dirichlet constraints for FEM");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_; }

    /*!
     * \brief Set the problem name.
     * \param newName The problem's name
     */
    void setName(const std::string& newName)
    {
        problemName_ = newName;
    }

    /*!
     * \name Boundary conditions and sources defining the problem
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be used for
     *        which equation on the degrees of freedom living on the given
     *        sub entity of the given grid element.
     * \param element The grid element on which the degree of freedom
     * \param subEntity The element's sub entity located on the domain boundary
     */
    template<class SubEntity>
    BoundaryTypes boundaryTypes(const Element& element, const SubEntity& subEntity) const
    { return asImp_().boundaryTypesAtPos(subEntity.geometry().center()); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation at a given position on the boundary.
     * \param globalPos The position on the boundary
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any boundaryTypes function
        //! set Dirichlet boundary conditions everywhere for all primary variables
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluate the Dirichlet boundary conditions for a boundary
     *        entity (sub entity of a grid element) that carries degrees of freedom.
     * \param element The grid element
     * \param subEntity The element's sub entity living on the boundary
     * \note The values will be assigned to all dofs living on that sub entity
     */
    template<class SubEntity>
    PrimaryVariables dirichlet(const Element& element, const SubEntity& subEntity) const
    { return asImp_().dirichletAtPos(subEntity.geometry().center()); }

    /*!
     * \brief Evaluate the Dirichlet boundary conditions at a given position on the boundary.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary segments are Dirichlet, "
                   "but does not provide a dirichlet() or a dirichletAtPos() function.");
    }

    /*!
     * \brief Enables / disables internal (non-boundary) Dirichlet constraints.
     * \todo TODO add internal Dirichlet constraint support for FEM
     */
    static constexpr bool enableInternalDirichletConstraints()
    { return false; }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     * \note This is the overload for the case where the Neumann condition is
     *       potentially solution dependent
     *
     * \param element The grid element
     * \param is The intersection
     * \param elemSol The element solution vector
     * \param ipData Shape function values/gradients evaluated at an integration point
     * \param ipVars The primary/secondary variables evaluated at an integration point
     */
    template<class ElementSolution, class IpData, class IpVariables>
    NumEqVector neumann(const Element& element,
                        const Intersection& is,
                        const ElementSolution& elemSol,
                        const IpData& ipData,
                        const IpVariables& ipVars) const
    { return asImp_().neumannAtPos(ipData.ipGlobal()); }

    /*!
     * \brief Evaluate the boundary conditions at a given position on a Neumann segment.
     * \param globalPos The position on a Neumann boundary segment
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any neumann function
        //! return no-flow Neumann boundary conditions at all Neumann boundaries
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the source term at a given integration point.
     * \note This is the overload for the case where the source term
     *       is potentially solution dependent and requires additional data.
     *
     * \param element The grid element
     * \param feGeometry The finite element geometry
     * \param elemSol The element solution vector
     * \param ipData Shape function values/gradients evaluated at an integration point
     * \param ipVars The primary/secondary variables evaluated at an integration point
     */
    template<class ElementSolution, class IpData, class IpVariables>
    NumEqVector source(const Element& element,
                       const FEElementGeometry& feGeometry,
                       const ElementSolution& elemSol,
                       const IpData& ipData,
                       const IpVariables& ipVars) const
    { return asImp_().sourceAtPos(ipData.ipGlobal()); }

    /*!
     * \brief Evaluate the source term at a given position.
     * \param globalPos The position of the integration point (in global coordinates)
     *                  at which the source term should be specified.
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any source function
        //! return 0.0 (no source terms)
        return NumEqVector(0.0);
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     * \param sol the initial solution vector
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

        sol.resize(gridGeometry_->numDofs());
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            const auto geometry = element.geometry();
            const auto refElement = ReferenceElements::general(geometry.type());

            auto feGeometry = localView(*gridGeometry_);
            feGeometry.bind(element);

            // loop over all dofs in this element and obtain initial condition
            const auto& fe = feGeometry.feBasisLocalView().tree().finiteElement();
            for (unsigned int i = 0; i < feGeometry.feBasisLocalView().size(); ++i)
            {
                const auto& localKey = fe.localCoefficients().localKey(i);
                auto subEntity = localKey.subEntity();
                auto codim = localKey.codim();

                // obtain the global position of the center of the sub-entity on
                // which the current local degree of freedom lives
                const auto globalPos = geometry.global(refElement.position(subEntity, codim));
                const auto dofIdxGlobal = feGeometry.feBasisLocalView().index(i);
                sol[dofIdxGlobal] = asImp_().initial(element, globalPos);
            }
        }
    }

    /*!
     * \brief Evaluate the initial conditions at the global position of a
     *        degree of freedom inside the given element.
     *
     * \param element The grid element
     * \param globalPos The position within the element
     */
    PrimaryVariables initial(const Element& element, const GlobalPosition& globalPos) const
    { return asImp_().initialAtPos(globalPos); }

    /*!
     * \brief Evaluate the initial conditions at a given position.
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no reasonable default value for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide an initial() or an initialAtPos() function.");
    }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     * \note This is the overload for extrusion being potentially solution-dependent.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     *
     * \param element The grid element
     * \param ipData Shape function values/gradients evaluated at an integration point
     * \param elemSol The element solution vector
     */
    template<class IpData, class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const IpData& ipData,
                           const ElementSolution& elemSol) const
    { return asImp_().extrusionFactorAtPos(ipData.ipGlobal()); }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        // As a default, i.e. if the user's problem does not overload
        // any extrusion factor function, return 1.0
        return 1.0;
    }

    // \}

    //! The grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The parameter group in which to retrieve runtime parameters
    const std::string& paramGroup() const
    { return paramGroup_; }

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    const Impl &asImp_() const { return *static_cast<const Impl *>(this); }

    //! \copydoc asImp_()
    Impl &asImp_() { return *static_cast<Impl *>(this); }

private:
    //! The grid geometry
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! The parameter group in which to retrieve runtime parameters
    std::string paramGroup_;

    //! The name of the problem
    std::string problemName_;
};

} // end namespace Dumux

#endif
