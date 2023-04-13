// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/grid/griddata.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>

namespace Dumux {

/*!
 * \brief Test problem for the (Navier-) Stokes model in a 3D channel
 *
 * Flow from left to right in a three-dimensional channel is considered. At the inlet (left)
 * and outlet (right) fixed values for pressure are set.
 */
template <class TypeTag, class BaseProblem>
class ThreeDChannelTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using DirichletValues = typename ParentType::DirichletValues;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    template<class GridData>
    ThreeDChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager, GridData gridData)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    , gridData_(gridData)
    {
        deltaP_ = getParam<Scalar>("Problem.DeltaP");
        density_ = getParam<Scalar>("Component.LiquidDensity");
        viscosity_ = getParam<Scalar>("Component.LiquidDynamicViscosity");
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (isOutlet_(scvf) || isInlet_(scvf))
                values.setAllNeumann();
            else
                values.setAllDirichlet();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scv The sub control volume
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolume& scv) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (fvGeometry.scv(scvf.insideScvIdx()).dofIndex() == scv.dofIndex() && scvf.boundary())
                {
                    if (isOutlet_(scvf) || isInlet_(scvf))
                        values.setAllNeumann();
                    else
                        values.setAllDirichlet();

                    break;
                }
            }
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // no-flow/no-slip
        return DirichletValues(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            if (isOutlet_(scvf) || isInlet_(scvf))
            {
                const auto p = isInlet_(scvf) ? deltaP_ : 0.0;
                values.axpy(-p, scvf.unitOuterNormal());
            }
        }
        else
        {
            if (isInlet_(scvf) || isOutlet_(scvf))
            {
                const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
                values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
            }
        }

        return values;
    }

    template<class ScvOrScvf>
    Scalar density(const Element&,
                   const FVElementGeometry&,
                   const ScvOrScvf&,
                   const bool isPreviousTimeStep = false) const
    {
        return density_;
    }

    template<class ScvOrScvf>
    Scalar effectiveViscosity(const Element&,
                              const FVElementGeometry&,
                              const ScvOrScvf&) const
    {
        return viscosity_;
    }

    template<class VelSolutionVector>
    void computeFluxes(const VelSolutionVector& sol)
    {
        Scalar influx = 0.0;
        Scalar outflux = 0.0;
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    if (isInlet_(scvf))
                        influx += scvf.area() * (sol[fvGeometry.scv(scvf.insideScvIdx()).dofIndex()] * scvf.unitOuterNormal());

                    else if (isOutlet_(scvf))
                        outflux += scvf.area() * (sol[fvGeometry.scv(scvf.insideScvIdx()).dofIndex()] * scvf.unitOuterNormal());
                }
            }
        }

        std::cout << "Influx: " << influx << " m^3/s" << std::endl;
        std::cout << "Outflux: " << outflux << " m^3/s" << std::endl;
        std::cout << "Balance: " << outflux+influx << " m^3/s" << std::endl;

        std::cout << "Transmissibility: " << outflux / deltaP_ << " m^3/(s*Pa)" << std::endl;
    }

#if HAVE_DUNE_FOAMGRID
    // Compute the wall shear stress
    template<class VelSolutionVector, class GridVariables, class CouplingManager, class Assembler>
    void computeWallShearStress(
        const VelSolutionVector& curSol,
        const GridVariables& gridVariables,
        CouplingManager& couplingManager,
        const Assembler& assembler
    ){
        const auto& gg = this->gridGeometry();
        const auto& gv = gg.gridView();

        using Grid = Dune::FoamGrid<dim-1, dimWorld>;
        Dune::GridFactory<Grid> factory;
        {
            std::vector<std::vector<std::vector<unsigned int>>> elems;
            elems.reserve(gg.numBoundaryScvf());
            std::vector<bool> insertedVertex(gv.size(dim), false);
            std::vector<int> vertexIndexMap(gv.size(dim), -1);
            std::vector<Dune::FieldVector<double, 3>> points;
            points.reserve(gg.numBoundaryScvf()*4);
            std::size_t boundaryElementIndex = 0;
            for (const auto& element : elements(gv))
            {
                for (const auto& intersection : intersections(gv, element))
                {
                    if (intersection.boundary())
                    {
                        const auto insideIdx = intersection.indexInInside();
                        const auto refElement = referenceElement(element);
                        const auto numVertices = refElement.size(insideIdx, 1, dim);
                        if (numVertices == 3)
                            elems.push_back({std::vector<unsigned int>{}});
                        else if (numVertices == 4)
                            elems.push_back({std::vector<unsigned int>{},
                                             std::vector<unsigned int>{}}); // add two triangles
                        else
                                DUNE_THROW(Dune::NotImplemented,
                                    "Wall shear stress for boundary type with " << numVertices << " corners"
                                );


                        for (int i = 0; i < numVertices; ++i)
                        {
                            const auto localVIdx = refElement.subEntity(insideIdx, 1, i, dim);
                            const auto& vertex = element.template subEntity<dim>(localVIdx);
                            const auto vIdx = gg.vertexMapper().index(vertex);
                            if (!insertedVertex[vIdx])
                            {
                                vertexIndexMap[vIdx] = points.size();
                                points.push_back(vertex.geometry().corner(0));
                                insertedVertex[vIdx] = true;
                            }


                            if (numVertices == 3)
                                elems[boundaryElementIndex][0].push_back(vertexIndexMap[vIdx]);

                            else if (numVertices == 4)
                            {
                                // add two triangles
                                if (i == 0)
                                    elems[boundaryElementIndex][0].push_back(vertexIndexMap[vIdx]);
                                else if (i == 1 || i == 2)
                                {
                                    elems[boundaryElementIndex][0].push_back(vertexIndexMap[vIdx]);
                                    elems[boundaryElementIndex][1].push_back(vertexIndexMap[vIdx]);
                                }
                                else
                                    elems[boundaryElementIndex][1].push_back(vertexIndexMap[vIdx]);
                            }
                            else
                                DUNE_THROW(Dune::NotImplemented,
                                    "Wall shear stress for boundary type with " << numVertices << " corners"
                                );
                        }

                        ++boundaryElementIndex;
                    }
                }
            }

            for (const auto& p : points)
                factory.insertVertex(p);
            for (const auto& e : elems)
                for (const auto& ee : e)
                    factory.insertElement(Dune::GeometryTypes::simplex(dim-1), ee);
        }

        auto bGrid = factory.createGrid();

        using ElementSet = GridViewGeometricEntitySet<typename Grid::LeafGridView, 0>;
        using Tree = BoundingBoxTree<ElementSet>;
        Tree tree(std::make_shared<ElementSet>(bGrid->leafGridView()));

        std::vector<GlobalPosition> wallShearStress(bGrid->leafGridView().size(0));

        auto fvGeometry = localView(gg);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        for (const auto& element : elements(gv))
        {
            const auto h = diameter(element.geometry());
            couplingManager.bindCouplingContext(
                CouplingManager::freeFlowMomentumIndex, element, assembler
            );
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    const auto ip = scvf.ipGlobal();

                    // assemble wss = sigma*n - (sigma*n*n)*n
                    const auto& fluxVarCache = elemFluxVarsCache[scvf];
                    using Tensor = Dune::FieldMatrix<Scalar, dim>;
                    Tensor sigma(0.0);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                       const auto& volVars = elemVolVars[scv];
                       for (int dir = 0; dir < dim; ++dir)
                           sigma[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(scv.indexInElement()));
                    }
                    sigma += getTransposed(sigma);
                    sigma *= -this->effectiveViscosity(element, fvGeometry, scvf);
                    for (int dir = 0; dir < dim; ++dir)
                        sigma[dir][dir] += this->pressure(element, fvGeometry, scvf);

                    const auto sigman = mv(sigma, scvf.unitOuterNormal());
                    const auto normalStress = (sigman * scvf.unitOuterNormal())*scvf.unitOuterNormal();
                    const auto wss = sigman - normalStress;

                    // make sure we hit all triangles of this face
                    const auto points = EmbeddedCoupling::circlePoints(ip, scvf.unitOuterNormal(), 0.05*h, 4);
                    for (const auto& p : points)
                        for (auto b : intersectingEntities(p, tree))
                            wallShearStress[b] = wss;
                }
            }
        }

        Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView> bMapper(
            bGrid->leafGridView(), Dune::mcmgLayout(Dune::Codim<0>())
        );
        using Field = Dumux::Vtk::Field<typename Grid::LeafGridView>;
        Dune::VTKWriter<typename Grid::LeafGridView> writer(bGrid->leafGridView());
        writer.addCellData(Field(
            bGrid->leafGridView(), bMapper, wallShearStress, "wallShearStress", dimWorld
        ).get());
        writer.write(this->name() + "_wall_shear_stress");
    }
#endif

private:
    bool isInlet_(const SubControlVolumeFace& scvf) const
    { return gridData_->getBoundaryDomainMarker(scvf.boundaryFlag()) == 1; }

    bool isOutlet_(const SubControlVolumeFace& scvf) const
    { return gridData_->getBoundaryDomainMarker(scvf.boundaryFlag()) == 2; }

    static constexpr Scalar eps_ = 1e-10;
    Scalar deltaP_, density_, viscosity_;

    std::shared_ptr<GridData<typename GridGeometry::GridView::Grid>> gridData_;
};

} // end namespace Dumux

#endif
