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
 *        Intersection computation relies on boundingBoxTree and therefore
          does not require grid-glue.
 */
#ifndef DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH
#define DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

//#include <dune/common/timer.hh>
//#include <dune/geometry/quadraturerules.hh>
//#include <dune/geometry/referenceelements.hh>
//#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/boundingboxtree.hh>
#include "couplingmapper.hh"
#include "couplingdata.hh"
#include <dumux/common/exceptions.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(StokesData);
NEW_PROP_TAG(DarcyData);
} // namespace Properties

/*!
 * \brief Manages the coupling between stokes elements and lower dimensional elements
 * \ingroup BoundaryCoupling
 */
template<class TypeTag>
class CouplingManagerStokesDarcy
{
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
        dimWorld = StokesGridView::dimensionworld
    };

    enum {
        stokesIsBox = GET_PROP_VALUE(StokesProblemTypeTag, ImplicitIsBox),
        darcyIsBox = GET_PROP_VALUE(DarcyProblemTypeTag, ImplicitIsBox)
    };

    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;

    using CoordScalar = typename StokesGridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using StokesVertex = typename StokesGridView::template Codim<stokesDim>::Entity;
    using DarcyVertex = typename DarcyGrid::template Codim<darcyDim>::Entity;

    using CouplingMapper = Dumux::CouplingMapperStokesDarcy<TypeTag>;
    using StokesData = typename GET_PROP_TYPE(TypeTag, StokesData);
    using DarcyData = typename GET_PROP_TYPE(TypeTag, DarcyData);

    using DofTypeIndices = typename GET_PROP(StokesProblemTypeTag, DofTypeIndices);
    using cellCenterIdx = typename DofTypeIndices::CellCenterIdx;
    using faceIdx = typename DofTypeIndices::FaceIdx;

public:

    /*!
     * \brief Constructor
     */
    CouplingManagerStokesDarcy(StokesProblem& stokesProblem, DarcyProblem& darcyProblem)
    : stokesProblem_(stokesProblem),
      darcyProblem_(darcyProblem),
      couplingMapper_(stokesProblem, darcyProblem, asImp_())
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
     * \brief Compute the coupling maps and stencil after the sub-problem have been initialized
     */
    void postInit()
    {
        couplingMapper_.computeCouplingMaps();
        computeStencils();
    }

    /*!
     * \brief Returns whether a Stokes element sub control volume is considered for coupling
     *        with the Darcy domain. This function is used for the boundaryType() method
     *        of the Stokes problem which handles vertices.
     *
     * \param element The element
     */
    bool isStokesCouplingEntity(const StokesElement &element) const
    {
        return !(couplingStencil(element).empty());
    }

    /*!
     * \brief Returns whether a Stokes element sub control volume face is considered for coupling
     *        with the Darcy domain. This function is used for the boundaryType() method
     *        of the Stokes problem which handles vertices.
     *
     * \param scvf The sub control volume face
     */
    bool isStokesCouplingEntity(const StokesSubControlVolumeFace& scvf) const
    {
        return couplingMapper_.stokesFaceToDarcyMap().count(scvf.dofIndex());
    }

    /*!
     * \brief Returns whether a Stokes element sub control volume is considered for coupling
     *        with the Darcy domain. This function is used for the boundaryType() method
     *        of the Stokes problem which handles vertices.
     *
     * \param scv The sub control volume
     */
    bool isStokesCouplingEntity(const StokesSubControlVolume& scv) const
    {
        return couplingMapper_.stokesCCToDarcyMap().count(scv.dofIndex());
    }

    /*!
     * \brief Returns whether a Darcy vertex is considered for coupling
     *        with the Stokes domain. This function is used for the boundaryType() method
     *        of the Darcy problem which handles vertices.
     *
     * \param vertex The vertex
     */
    bool isDarcyCouplingEntity(const DarcyVertex &vertex) const
    {
        const auto darcyDofIdxGlobal = darcyGridView().indexSet().index(vertex);
        return couplingMapper_.darcyToStokesMap().count(darcyDofIdxGlobal);
    }

    /*!
     * \brief Returns whether a Darcy vertex is considered for coupling
     *        with the Stokes domain. This function is used for the boundaryType() method
     *        of the Darcy problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcyElement &element) const
    {
        return !(couplingStencil(element).empty());
    }

    /*!
     * \brief Returns whether a Darcy vertex is considered for coupling
     *        with the Stokes domain. This function is used for the boundaryType() method
     *        of the Darcy problem which handles vertices.
     *
     * \param element The element
     */
    bool isDarcyCouplingEntity(const DarcySubControlVolume &scv) const
    {
        return couplingMapper_.darcyToStokesMap().count(scv.dofIndex());
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const CouplingMapper &couplingMapper() const
    {
        return couplingMapper_;
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &darcyToStokesMap() const
    {
        return couplingMapper_.darcyToStokesMap();
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &stokesFaceToDarcyMap() const
    {
        return couplingMapper_.stokesFaceToDarcyMap();
    }

    /*!
     * \brief Returns a reference to the coupling mapper.
     */
    const auto &stokesCCToDarcyMap() const
    {
        return couplingMapper_.stokesCCToDarcyMap();
    }


    void computeStencils()
    {
        const auto &stokesTree = stokesProblem_.boundingBoxTree();

        // compute Stokes cell-center coupling stencil using the coupling map
        for (const auto &entry : couplingMapper_.stokesCCToDarcyMap())
        {
            const auto stokesElementIdx = entry.first;
            Dune::dverb << "Stokes element " << stokesElementIdx <<  " is coupled to Darcy element:";

            // get the Darcy DOFs associated with the Stokes element
            auto &darcyInfo = couplingMapper_.stokesCCToDarcyMap().at(stokesElementIdx);
            auto darcyElementIndex = darcyInfo.darcyElementIdx;
            Dune::dverb << " " << darcyElementIndex ;
            stokesCCCouplingStencils_[stokesElementIdx].push_back(darcyElementIndex);

            const auto& stokesElement = stokesTree.entity(stokesElementIdx);
            StokesFVElementGeometry stokesFvGeometry = localView(stokesProblem_.model().globalFvGeometry());
            stokesFvGeometry.bind(stokesElement);

            for(auto&& stokesScvf : scvfs(stokesFvGeometry))
                stokesFaceCouplingStencils_[stokesScvf.dofIndex()].push_back(darcyElementIndex);

            Dune::dverb << std::endl;
        }

        Dune::dverb << std::endl;

        // compute Darcy coupling stencil using the coupling map
        for (const auto &entry : couplingMapper_.darcyToStokesMap())
        {
            const auto darcyElementIdx = entry.first;
            Dune::dverb << "Darcy element " << darcyElementIdx <<  " is coupled to Stokes element:";

            // get the Stokes DOFs associated with the Darcy element
            const auto &stokesElements = couplingMapper_.darcyToStokesMap().at(darcyElementIdx).stokesElementIndices();
            auto stokesElementIndex = stokesElements[0];
            Dune::dverb << " Stokes element " << stokesElementIndex;
            darcyToCCCouplingStencils_[darcyElementIdx].push_back(stokesElementIndex);

            const auto& stokesElement = stokesTree.entity(stokesElementIndex);
            StokesFVElementGeometry stokesFvGeometry = localView(stokesProblem_.model().globalFvGeometry());
            stokesFvGeometry.bind(stokesElement);

            for(auto&& stokesScvf : scvfs(stokesFvGeometry))
                darcyToFaceCouplingStencils_[darcyElementIdx].push_back(stokesScvf.dofIndex());

            Dune::dverb << std::endl;
        }

        auto sortAndMakeUnique = [](auto&& stencils)
        {
            for(auto&& stencil : stencils)
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        };

        // sort and make unique
        sortAndMakeUnique(stokesCCCouplingStencils_);
        sortAndMakeUnique(stokesFaceCouplingStencils_);
        sortAndMakeUnique(darcyToCCCouplingStencils_);
        sortAndMakeUnique(darcyToFaceCouplingStencils_);
    }

    /*!
     * \brief Returns a coupling stencil for a Stokes element
     */
    const std::vector<unsigned int>& couplingStencil(const StokesElement& element) const
    {
        const unsigned int eIdx = stokesProblem().elementMapper().index(element);
        if (stokesCCCouplingStencils_.count(eIdx))
            return stokesCCCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Stokes element
     */
    const std::vector<unsigned int>& couplingStencil(const StokesSubControlVolumeFace& scvf) const
    {
        const unsigned int dofIdx = scvf.dofIndex();
        if (stokesFaceCouplingStencils_.count(dofIdx))
            return stokesFaceCouplingStencils_.at(dofIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Darcy element
     */
    const std::vector<unsigned int>& couplingStencil(const DarcyElement& element, cellCenterIdx) const
    {
        const unsigned int eIdx = darcyProblem().elementMapper().index(element);
        if (darcyToCCCouplingStencils_.count(eIdx))
            return darcyToCCCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a coupling stencil for a Darcy element
     */
    const std::vector<unsigned int>& couplingStencil(const DarcyElement& element, faceIdx) const
    {
        const unsigned int eIdx = darcyProblem().elementMapper().index(element);
        if (darcyToFaceCouplingStencils_.count(eIdx))
            return darcyToFaceCouplingStencils_.at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief Returns a class containing all relevant Stokes data
     */
    auto stokesData() const
    {
        StokesData data(*this);
        return data;
    }

    /*!
     * \brief Returns a class containing all relevant Darcy data
     */
    auto darcyData() const
    {
        DarcyData data(*this);
        return data;
    }

    //! evaluate coupling residual for the derivative stokes DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
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
        // calculate the local residual of the Stokes element
        auto prevElemVolVars = localView(stokesProblem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, stokesProblem().model().prevSol());
        stokesLocalResidual_.evalCellCenter(element, fvGeometry, scvI,
                                          prevElemVolVars, curElemVolVars,
                                          prevGlobalFaceVars, curGlobalFaceVars,
                                          elemBcTypes, elemFluxVarsCache);
        return stokesLocalResidual_.ccResidual();
    }

    //! evaluate coupling residual for the derivative stokes DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    auto evalStokesFaceCouplingResidual(const StokesElement& element,
                              const StokesFVElementGeometry& fvGeometry,
                              const StokesSubControlVolumeFace& scvf,
                              const StokesElementVolumeVariables& curElemVolVars,
                              const StokesGlobalFaceVars& prevGlobalFaceVars,
                              StokesGlobalFaceVars& curGlobalFaceVars,
                              const StokesElementBoundaryTypes& elemBcTypes,
                              const StokesElementFluxVariablesCache& elemFluxVarsCache)
    {
        // calculate the local residual of the Stokes element
        auto prevElemVolVars = localView(stokesProblem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, stokesProblem().model().prevSol());
        stokesLocalResidual_.evalFace(element, fvGeometry, scvf,
                                    prevElemVolVars, curElemVolVars,
                                    prevGlobalFaceVars, curGlobalFaceVars,
                                    elemBcTypes, elemFluxVarsCache, true /*resizeResidual*/);
        return stokesLocalResidual_.faceResidual(scvf.localFaceIdx());
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to stokes DOF
    //! we only need to evaluate the part of the residual that will be influence by the stokes DOF
    auto evalDarcyCouplingResidual(const DarcyElement& element,
                              const DarcyFVElementGeometry& fvGeometry,
                              const DarcyElementVolumeVariables& curElemVolVars,
                              const DarcyElementBoundaryTypes& elemBcTypes,
                              const DarcyElementFluxVariablesCache& elemFluxVarsCache)
    {
        // calculate the local residual of the Darcy element
        auto prevElemVolVars = localView(darcyProblem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, darcyProblem().model().prevSol());
        darcyLocalResidual_.eval(element, fvGeometry, prevElemVolVars, curElemVolVars, elemBcTypes, elemFluxVarsCache);
        return darcyLocalResidual_.residual();
    }

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
    CouplingMapper couplingMapper_;

    mutable StokesLocalResidual stokesLocalResidual_;
    mutable DarcyLocalResidual darcyLocalResidual_;

    std::unordered_map<unsigned int, std::vector<unsigned int> > stokesCCCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > stokesFaceCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToCCCouplingStencils_;
    std::unordered_map<unsigned int, std::vector<unsigned int> > darcyToFaceCouplingStencils_;
    std::vector<unsigned int> emptyStencil_;
};

} // namespace Dumux

#endif // DUMUX_COUPLINGMANAGER_STOKES_DARCY_HH
