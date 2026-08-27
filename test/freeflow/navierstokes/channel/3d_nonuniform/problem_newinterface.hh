// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_NEWINTERFACE_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_NEWINTERFACE_HH

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/grid/griddata.hh>
#include <dumux/common/boundaryflag.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/multidomain/embedded/circlepoints.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/typetraits/griddiscretization.hh>

namespace Dumux {

/*!
 * \brief Test problem for the (Navier-) Stokes model in a 3D channel
 *
 * Flow from left to right in a three-dimensional channel is considered. At the inlet (left)
 * and outlet (right) fixed values for pressure are set.
 */
template <class TypeTag, class BaseProblem>
class ThreeDChannelTestProblemNewInterface : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using SubControlVolume = typename ElementDiscretization::SubControlVolume;
    using SubControlVolumeFace = typename ElementDiscretization::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using BoundaryFlag = Dumux::BoundaryFlag<typename GridDiscretization::GridView::Grid>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using DirichletValues = typename ParentType::DirichletValues;
    using GridIndexType = typename IndexTraits<typename GridDiscretization::GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int dim = GridDiscretization::GridView::dimension;
    static constexpr auto dimWorld = GridDiscretization::GridView::dimensionworld;

    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    template<class GridData>
    ThreeDChannelTestProblemNewInterface(std::shared_ptr<const GridDiscretization> gridDiscretization, std::shared_ptr<CouplingManager> couplingManager, GridData gridData)
    : ParentType(gridDiscretization, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    , gridData_(gridData)
    {
        deltaP_ = getParam<Scalar>("Problem.DeltaP");
        density_ = getParam<Scalar>("Component.LiquidDensity");
        viscosity_ = getParam<Scalar>("Component.LiquidDynamicViscosity");

        if constexpr (ParentType::isMomentumProblem())
            appendDirichletConstraints_();
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary face.
     *
     * \param elemDisc The element discretization
     * \param boundaryFace The boundary face
     */
    BoundaryTypes boundaryTypes(const ElementDiscretization& elemDisc,
                                const ElementDiscretization::BoundaryFace& boundaryFace) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto flag = boundaryFace.boundaryFlag();
            if (isMomentumFluxBoundary_(flag))
                values.setAllFluxBoundary();
        }
        else
            values.setFluxBoundary(Indices::conti0EqIdx);

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

    const auto& constraints() const
    { return constraints_; }

    /*!
     * \brief Evaluates the boundary flux related to a localDof at a given interpolation point.
     *
     * \param elemDisc The element discretization
     * \param elemVars All variables related to the element
     * \param faceIpData Interpolation point data
     */
    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    {
        BoundaryFluxes values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            if (isMomentumFluxBoundary_(boundaryFlag_(elemDisc, faceIpData)))
            {
                const auto p = isInlet_(boundaryFlag_(elemDisc, faceIpData)) ? deltaP_ : 0.0;
                values.axpy(-p, faceIpData.unitOuterNormal());
            }
        }
        else
        {
            const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
            if (isInlet_(scvf.boundaryFlag()) || isOutlet_(scvf.boundaryFlag()))
            {
                const auto insideDensity = elemVars[scvf.insideScvIdx()].density();
                values[Indices::conti0EqIdx] = this->velocity(elemDisc, faceIpData) * insideDensity * scvf.unitOuterNormal();
            }
        }

        return values;
    }

    template<class ScvOrScvfOrIpData>
    Scalar density(const Element&,
                   const ElementDiscretization&,
                   const ScvOrScvfOrIpData&,
                   const bool isPreviousTimeStep = false) const
    {
        return density_;
    }

    template<class ScvOrScvfOrIpData>
    Scalar effectiveViscosity(const Element&,
                              const ElementDiscretization&,
                              const ScvOrScvfOrIpData&) const
    {
        return viscosity_;
    }

    template<class VelSolutionVector>
    void computeFluxes(const VelSolutionVector& sol)
    {
        Scalar influx = 0.0;
        Scalar outflux = 0.0;
        const auto& gridDiscretization = Dumux::gridDiscretization(*this);
        auto elemDisc = localView(gridDiscretization);
        for (const auto& element : elements(gridDiscretization.gridView()))
        {
            elemDisc.bind(element);
            for (const auto& scvf : scvfs(elemDisc))
            {
                if (scvf.boundary())
                {
                    if (isInlet_(scvf.boundaryFlag()))
                        influx += scvf.area() * (sol[elemDisc.scv(scvf.insideScvIdx()).dofIndex()] * scvf.unitOuterNormal());

                    else if (isOutlet_(scvf.boundaryFlag()))
                        outflux += scvf.area() * (sol[elemDisc.scv(scvf.insideScvIdx()).dofIndex()] * scvf.unitOuterNormal());
                }
            }
        }

        std::cout << "Influx: " << influx << " m^3/s" << std::endl;
        std::cout << "Outflux: " << outflux << " m^3/s" << std::endl;
        std::cout << "Balance: " << outflux+influx << " m^3/s" << std::endl;

        std::cout << "Transmissibility: " << outflux / deltaP_ << " m^3/(s*Pa)" << std::endl;
    }

private:
    template<class IpData>
    auto boundaryFlag_(const ElementDiscretization& elemDisc, const IpData& ipData) const
    {
        if constexpr (Dumux::Concept::ScvfIpData<IpData>)
            return elemDisc.scvf(ipData.scvfIndex()).boundaryFlag();
        else if constexpr (Dumux::Concept::BoundaryFaceIpData<IpData>)
            return elemDisc.boundaryFace(ipData.boundaryFaceIndex()).boundaryFlag();
        else
            DUNE_THROW(Dune::InvalidStateException, "Unknown type of interpolation point data!");
    }

    template<class BoundaryFlag>
    bool isMomentumFluxBoundary_(const BoundaryFlag& bFlag) const
    { return (isOutlet_(bFlag) || isInlet_(bFlag)); }

    void appendDirichletConstraints_()
    {
        const auto& gridDiscretization = Dumux::gridDiscretization(*this);
        auto elemDisc = localView(gridDiscretization);
        for (const auto& element : elements(gridDiscretization.gridView()))
        {
            elemDisc.bind(element);

            for(const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if(!isMomentumFluxBoundary_(boundaryFace.boundaryFlag()))
                {
                    for(const auto& localDof : localDofs(elemDisc, boundaryFace))
                    {
                        const auto globalPos = ipData(elemDisc, localDof).global();
                        ConstraintInfo info;
                        info.setAll();

                        DirichletValues dirichletValues(this->dirichletAtPos(globalPos));
                        constraints_.push_back(DirichletConstraintData{std::move(info), std::move(dirichletValues), localDof.dofIndex()});
                    }
                }
            }
        }
    }

    bool isInlet_(int flag) const
    { return gridData_->getBoundaryDomainMarker(flag) == 1; }

    bool isOutlet_(int flag) const
    { return gridData_->getBoundaryDomainMarker(flag) == 2; }

    static constexpr Scalar eps_ = 1e-10;
    Scalar deltaP_, density_, viscosity_;

    std::shared_ptr<GridData<typename GridDiscretization::GridView::Grid>> gridData_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
