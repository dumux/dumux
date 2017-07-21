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
 * \ingroup BoundaryCoupling
 * \brief Coupling manager for two domains with the same dimension.
 */
#ifndef DUMUX_COUPLINGMANAGER_HH
#define DUMUX_COUPLINGMANAGER_HH

// TODO clean up
#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

// #include <dumux/common/boundingboxtree.hh>
// #include "couplingmapper.hh"
// #include "couplingdata.hh"
#include <dumux/common/exceptions.hh>

namespace Dumux
{

namespace Properties // TODO clean up
{
// Property forward declarations
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(StokesData);
NEW_PROP_TAG(DarcyData);
} // namespace Properties

/*!
 * \brief Manages the coupling between staggered grid and cell-centered finite volume elements
 * \ingroup BoundaryCoupling
 */
template<class TypeTag>
class BoundaryCouplingManagerPNMStokes
{
    // TODO clean up typedefs
    using Implementation = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    using StokesGrid = typename GET_PROP_TYPE(StokesProblemTypeTag, Grid);
    using DarcyGrid = typename GET_PROP_TYPE(DarcyProblemTypeTag, Grid);

    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using DarcySubControlVolume = typename GET_PROP_TYPE(DarcyProblemTypeTag, SubControlVolume);
    using StokesSubControlVolume = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolume);
    using StokesSubControlVolumeFace = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolumeFace);

    using DarcyPrimaryVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, PrimaryVariables);
    using StokesPrimaryVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, PrimaryVariables);

    using DarcyElementVolumeVariables = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementVolumeVariables);
    using StokesElementVolumeVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementVolumeVariables);

    using StokesElementBoundaryTypes = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementBoundaryTypes);
    using StokesElementFluxVariablesCache = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementFluxVariablesCache);

    using DarcyElementBoundaryTypes = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementBoundaryTypes);
    using DarcyElementFluxVariablesCache = typename GET_PROP_TYPE(DarcyProblemTypeTag, ElementFluxVariablesCache);

    using DarcyFVElementGeometry = typename GET_PROP_TYPE(DarcyProblemTypeTag, FVElementGeometry);
    using StokesFVElementGeometry = typename GET_PROP_TYPE(StokesProblemTypeTag, FVElementGeometry);

    using StokesLocalResidual = typename GET_PROP_TYPE(StokesProblemTypeTag, LocalResidual);
    using DarcyLocalResidual = typename GET_PROP_TYPE(DarcyProblemTypeTag, LocalResidual);

    using StokesGlobalFaceVars = typename GET_PROP_TYPE(StokesProblemTypeTag, GlobalFaceVars);

    enum {
        stokesDim = StokesGridView::dimension,
        darcyDim = DarcyGridView::dimension,
//         dim = StokesGridView::dimension,
        dimWorld = StokesGridView::dimensionworld
    };

    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;

    using CoordScalar = typename StokesGridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using StokesVertex = typename StokesGridView::template Codim<stokesDim>::Entity; // TODO
    using DarcyVertex = typename DarcyGrid::template Codim<darcyDim>::Entity; // TODO

    using CouplingMapper = Dumux::BoundaryCouplingMapperPNMStokes<TypeTag>; // TODO
    using StokesData = typename GET_PROP_TYPE(TypeTag, StokesData);
    using DarcyData = typename GET_PROP_TYPE(TypeTag, DarcyData);

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    using cellCenterIdx = typename DofTypeIndices::CellCenterIdx;
    using faceIdx = typename DofTypeIndices::FaceIdx;

public:

    // TODO move duplicate code into function (determine coupling interface, ...)


    /*!
     * \brief Constructor
     */
    CouplingManagerStokesDarcy(StokesProblem& stokesProblem, DarcyProblem& darcyProblem)
    : stokesProblem_(stokesProblem),
      darcyProblem_(darcyProblem),
      couplingMapper_(stokesProblem, darcyProblem, asImp_()) // TODO without couplingmapper!
    {
        // initialize the local residuals
        stokesLocalResidual_.init(stokesProblem);
        darcyLocalResidual_.init(darcyProblem);
    }

    /*!
     * \brief Return a reference to the Stokes problem
     */
    const StokesProblem& stokesProblem() const
    {
        return stokesProblem_;
    }

    /*!
     * \brief Return a reference to the Darcy problem
     */
    const DarcyProblem& darcyProblem() const
    {
        return darcyProblem_;
    }

    /*!
     * \brief Return a reference to the Stokes gridview
     */
    const StokesGridView& stokesGridView() const
    {
        return stokesProblem().gridView();
    }

    /*!
     * \brief Return a reference to the Darcy gridview
     */
    const DarcyGridView& darcyGridView() const
    {
        return darcyProblem().gridView();
    }

    void preInit()
    {}


    /*!
     * \brief Compute the stencils after the subproblems have been initialized
     */
    void postInit()
    {
//         couplingMapper_.computeCouplingMaps(); // TODO?
        computeStencils(); // TODO ?
    }



    void computeStencils() // TODO?? --> see staggered grid!
    {

        // compute stokes/stokes cell-center coupling stencil using the coupling map


        // compute lowdim/darcy coupling stencil using the coupling map


        // sort and make unique
    }

     /*!
     * \brief Returns a class containing all relevant Stokes data
     */
    auto stokesData() const // TODO necessary?
    {
        StokesData data(*this);
        return data;
    }

    /*!
     * \brief Returns a class containing all relevant Darcy data
     */
    auto darcyData() const // TODO necessary?
    {
        DarcyData data(*this);
        return data;
    }


    //! evaluate coupling residual for the derivative stokes DOF with respect to darcy DOF
    // TODO -- different structure? (q^pm = q^ff?)
    auto evalStokesCCCouplingResidual(const StokesElement& element,
                                      const StokesFVElementGeometry& fvGeometry,
                                      const StokesElementVolumeVariables& curElemVolVars,
                                      const StokesGlobalFaceVars& prevGlobalFaceVars,
                                      StokesGlobalFaceVars& curGlobalFaceVars,
                                      const StokesElementBoundaryTypes& elemBcTypes,
                                      const StokesElementFluxVariablesCache& elemFluxVarsCache)
    {
        auto ccGlobalI =  stokesProblem().elementMapper().index(element);
        auto&& scvI = fvGeometry.scv(ccGlobalI);

        // calculate the local residual of the stokes element
        auto prevElemVolVars = localView(stokesProblem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, stokesProblem().model().prevSol());
//         stokesLocalResidual_.evalCellCenter(element, fvGeometry, scvI,
//                                           prevElemVolVars, curElemVolVars,
//                                           prevGlobalFaceVars, curGlobalFaceVars,
//                                           elemBcTypes, elemFluxVarsCache);
//      return stokesLocalResidual_.ccResidual();

        // assumption: method only called for elements which are at the interface
//         for (auto scvf : fvGeometry.scvfs(fvGeometry)) // TODO ??? -- rough draft
//         {
//             if (abs(scvf.center()[1] - problem.interfaceVerticalPos()) < eps_) // |y(scvf.center) - y(interface)| < eps ???
//             {
//                 auto couplingScfvIdx = scvf.index();
//                 continue;
//             }
//         }

        auto outerUnitNormal = fvGeometry.scfv(couplingScfvIdx).unitOuterNormal();
        Scalar interfaceArea = fvGeometry.scfv(couplingScfvIdx).area();

        Scalar densityStokes = prevElemVolVars.density();
        const Scalar velocity = prevGlobalFaceVars.faceVars(scvf.dofIndex()).velocity();
        auto normalVelocityAtInterface = velocity * outerUnitNormal;

        // mass coupling condition
        auto stokesCCCouplingResidual = densityStokes * normalVelocityAtInterface * interfaceArea; // TODO nc --> rho_g^ff

        return stokesCCCouplingResidual;
    }


    auto evalStokesFaceCouplingResidual(const StokesElement& element,
                                        const StokesFVElementGeometry& fvGeometry,
                                        const StokesSubControlVolumeFace& scvf,
                                        const StokesElementVolumeVariables& curElemVolVars,
                                        const StokesGlobalFaceVars& prevGlobalFaceVars,
                                        StokesGlobalFaceVars& curGlobalFaceVars,
                                        const StokesElementBoundaryTypes& elemBcTypes,
                                        const StokesElementFluxVariablesCache& elemFluxVarsCache)
    {
        // calculate the local residual of the stokes element
        auto prevStokesElemVolVars = localView(stokesProblem().model().prevGlobalVolVars());
        prevStokesElemVolVars.bindElement(element, fvGeometry, stokesProblem().model().prevSol());
        auto prevDarcyElemVolVars = localView(darcyProblem().model().prevGlobalVolVars());

//         stokesLocalResidual_.evalFace(element, fvGeometry, scvf,
//                                     prevElemVolVars, curElemVolVars,
//                                     prevGlobalFaceVars, curGlobalFaceVars,
//                                     elemBcTypes, elemFluxVarsCache, true /*resizeResidual*/);
//         return stokesLocalResidual_.faceResidual(scvf.localFaceIdx());

        // assumption: method only called for elements which are at the interface
//         for (auto scvf : fvGeometry.scvfs(fvGeometry)) // TODO ??? -- rough draft
//         {
//             if (abs(scvf.center()[1] - problem.interfaceVerticalPos()) < eps_)
//             {
//                 auto couplingScfvIdx = scvf.index();
//                 continue;
//             }
//         }

        auto outerUnitNormal = fvGeometry.scfv(couplingScfvIdx).unitOuterNormal();
        Scalar interfaceArea = fvGeometry.scfv(couplingScfvIdx).area();
        Scalar pressureStokes = prevDarcyElemVolVars.pressure();

        // Scalar permeability = spatialParams_.intrinsicPermeabilityAtPos(faceCenterGlobal);
        Scalar sqrtPerm = std::sqrt(prevDarcyElemVolVars.permeability());
        Scalar alphaBJ = //  Scalar alphaBeaversJoseph = spatialParams_.beaversJosephCoeffAtPos(faceCenterGlobal);
        Scalar centerDistance = // elementCenetersGlobal - faceCenterGlobal
        Scalar beta = -1.0 * sqrtPerm / alphaBJ / centerDistance * unitOuterNormal;

        // normal momentum coupling condition ! originally formulated for pm ! TODO --> pressureDarcy, * -1 ???
        auto stokesFaceCouplingResidual = pressureStokes * outerUnitNormal * interfaceArea;

        // tangential momentoum coupling condition --> Dirichlet b.c.??
        std::vector<int> tangDims {0, 1}; // TODO 2D !
        for (auto curTangDim : tangDims)
        {
            Scalar BJvelocity0 = beta * velocity[curTangDim*2][curTangDim] / (1.0 + beta);
            Scalar BJvelocity1 = beta * velocity[curTangDim*2 + 1][curTangDim] / (1.0 + beta);

            Scalar tangMomResidual =

            stokesFaceCouplingResidual += tangMomResidual;
        }

        return stokesFaceCouplingResidual;
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to stokes DOF
    // TODO

    // auto evalCouplingResidual(const DarcyElement& element, ...






    // q^if berechnen (mit ff oder pm) --> wo addieren? StokesCC, StokesFace, DarcyCC??






    protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    StokesProblem& stokesProblem_;
    DarcyProblem& darcyProblem_;
//     CouplingMapper couplingMapper_;

    mutable StokesLocalResidual stokesLocalResidual_;
    mutable DarcyLocalResidual darcyLocalResidual_;

    Scalar eps_ = 1e-9; // TODO ?

//     std::unordered_map<unsigned int, std::vector<unsigned int> > stokesCCCouplingStencils_;
//     std::unordered_map<unsigned int, std::vector<unsigned int> > stokesFaceCouplingStencils_;
//     std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToCCCouplingStencils_;
//     std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToFaceCouplingStencils_;
//     std::vector<unsigned int> emptyStencil_;

};

} // namespace Dumux

#endif