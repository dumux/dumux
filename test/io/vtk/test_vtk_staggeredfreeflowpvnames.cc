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
 * \brief Test for writing and reading vtk files for the staggered grid free flow models
 */
#include <config.h>

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/compositional/kepsilonncmodel.hh>
#include <dumux/freeflow/compositional/lowrekepsilonncmodel.hh>
#include <dumux/freeflow/compositional/komegancmodel.hh>
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/compositional/zeroeqncmodel.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/loadsolution.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>

namespace Dumux {
namespace Properties {

NEW_TYPE_TAG(StaggeredPVNamesTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel));

NEW_TYPE_TAG(NavierStokesPVNameTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, NavierStokes));
NEW_TYPE_TAG(NavierStokesNIPVNameTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, NavierStokesNI));
NEW_TYPE_TAG(NavierStokesNCPVNameTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, NavierStokesNC));
NEW_TYPE_TAG(NavierStokesNCNIPVNameTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, NavierStokesNCNI));

NEW_TYPE_TAG(KEpsilonNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KEpsilon));
NEW_TYPE_TAG(KEpsilonNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KEpsilonNI));
NEW_TYPE_TAG(KEpsilonNCNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KEpsilonNC));
NEW_TYPE_TAG(KEpsilonNCNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KEpsilonNCNI));

NEW_TYPE_TAG(LowReKEpsilonNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, LowReKEpsilon));
NEW_TYPE_TAG(LowReKEpsilonNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, LowReKEpsilonNI));
NEW_TYPE_TAG(LowReKEpsilonNCNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, LowReKEpsilonNC));
NEW_TYPE_TAG(LowReKEpsilonNCNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, LowReKEpsilonNCNI));

NEW_TYPE_TAG(KOmegaNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KOmega));
NEW_TYPE_TAG(KOmegaNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KOmegaNI));
NEW_TYPE_TAG(KOmegaNCNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KOmegaNC));
NEW_TYPE_TAG(KOmegaNCNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, KOmegaNCNI));

NEW_TYPE_TAG(OneEqNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, OneEq));
NEW_TYPE_TAG(OneEqNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, OneEqNI));
NEW_TYPE_TAG(OneEqNCNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, OneEqNC));
NEW_TYPE_TAG(OneEqNCNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, OneEqNCNI));

NEW_TYPE_TAG(ZeroEqNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, ZeroEq));
NEW_TYPE_TAG(ZeroEqNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, ZeroEqNI));
NEW_TYPE_TAG(ZeroEqNCNameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, ZeroEqNC));
NEW_TYPE_TAG(ZeroEqNCNINameTestTypeTag, INHERITS_FROM(StaggeredPVNamesTestTypeTag, ZeroEqNCNI));

// The fluid system
SET_PROP(StaggeredPVNamesTestTypeTag, FluidSystem)
{
  using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

SET_TYPE_PROP(StaggeredPVNamesTestTypeTag, Scalar, double);

// Set the grid type
SET_TYPE_PROP(StaggeredPVNamesTestTypeTag, Grid, Dune::YaspGrid<2>);

SET_PROP(StaggeredPVNamesTestTypeTag, Problem)
{
private:
    // use the ZeroEqProblem as base class for non-RANS models and for the ZeroEq model
    // use the the KEpsilonProblem as base class for all RANS models except the ZeroEq model
    // NOTE: this rather unpleasant hack will be removed once the RANS models have been unified
    using MTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static constexpr auto dim = MTraits::dim();
    static constexpr auto nComp = MTraits::numComponents();
    static constexpr auto numEq = MTraits::numEq();

    using BaseTurbulentProblem = std::conditional_t<(std::is_same<typename MTraits::Indices, KOmegaIndices<dim, nComp>>::value ||
                                                     std::is_same<typename MTraits::Indices, FreeflowNonIsothermalIndices<KOmegaIndices<dim, nComp>, numEq>>::value),
                                                     KOmegaProblem<TypeTag>,
                                                     std::conditional_t<(std::is_same<typename MTraits::Indices, KEpsilonIndices<dim, nComp>>::value ||
                                                                         std::is_same<typename MTraits::Indices, FreeflowNonIsothermalIndices<KEpsilonIndices<dim, nComp>, numEq>>::value),
                                                                         KEpsilonProblem<TypeTag>,
                                                                         std::conditional_t<(std::is_same<typename MTraits::Indices, LowReKEpsilonIndices<dim, nComp>>::value ||
                                                                                             std::is_same<typename MTraits::Indices, FreeflowNonIsothermalIndices<LowReKEpsilonIndices<dim, nComp>, numEq>>::value),
                                                                         LowReKEpsilonProblem<TypeTag>, OneEqProblem<TypeTag>>>>;

    using BaseProblem = std::conditional_t<MTraits::usesTurbulenceModel() &&
                                           !std::is_same<MTraits, RANSModelTraits<dim>>::value &&
                                           !std::is_same<MTraits, FreeflowNIModelTraits<RANSModelTraits<dim>>>::value &&
                                           !std::is_same<MTraits, ZeroEqNCModelTraits<dim, nComp, false, 0>>::value &&
                                           !std::is_same<MTraits, FreeflowNIModelTraits<ZeroEqNCModelTraits<dim, nComp, false, 0>>>::value,
                                           BaseTurbulentProblem, ZeroEqProblem<TypeTag>>;

    template<class TTag>
    class MockProblem : public BaseProblem
    {
        using ParentType = BaseProblem;
        using BoundaryTypes = typename GET_PROP_TYPE(TTag, BoundaryTypes);
        using Scalar = typename GET_PROP_TYPE(TTag, Scalar);
        using Traits = typename GET_PROP_TYPE(TTag, ModelTraits);
    public:
        using ParentType::ParentType;

        template<class GlobalPosition>
        BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
        {
             BoundaryTypes values;
             return values;
        }

        Scalar temperature() const
        { return 300; }

        template<class T = TTag, bool enable = GET_PROP_TYPE(T, ModelTraits)::usesTurbulenceModel(), std::enable_if_t<!enable, int> = 0>
        void updateStaticWallProperties() {}

        template<class U, bool enable = Traits::usesTurbulenceModel(), std::enable_if_t<!enable, int> = 0>
        void updateDynamicWallProperties(const U&) {}

        // for ZeroEq model
        template<class T = TTag, bool enable = GET_PROP_TYPE(T, ModelTraits)::usesTurbulenceModel(), std::enable_if_t<enable, int> = 0>
        void updateStaticWallProperties()
        { ParentType::updateStaticWallProperties(); }

        // for ZeroEq model
        template<class U, bool enable = Traits::usesTurbulenceModel(), std::enable_if_t<enable, int> = 0>
        void updateDynamicWallProperties(const U& u)
        { return ParentType::updateDynamicWallProperties(u); }

        template<class Scvf>
        bool isOnWall(const Scvf&) const
        { return true; }
    };

public:
    using type = MockProblem<TypeTag>;
};

}
}

template<class SolutionVector, class Values>
void assignValues(SolutionVector& sol, Values values)
{
    for (auto& entry : sol)
    {
        for (int pvIdx = 0; pvIdx < values.size(); ++pvIdx)
        {
            // make sure to get values that can be exactly represented in the vtk file (Float32)
            std::stringstream stream;
            stream << std::fixed << std::setprecision(5) << std::scientific << values[pvIdx];
            entry[pvIdx] = std::stod(stream.str());
        }

        // increment all values by 1% for each dof
        for (auto& i : values)
            i += 0.01*i;
    }
}

template<class TypeTag, class FVGridGeometry, std::size_t numValues>
void testWriteAndReadVtk(std::shared_ptr<FVGridGeometry> fvGridGeometry,
                         const std::array<typename GET_PROP_TYPE(TypeTag, Scalar), numValues>& values,
                         const std::string& fileName,
                         bool verbose = false,
                         bool deleteFiles = true)
{
    using namespace Dumux;

    // the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    SolutionVector writeFrom;

    writeFrom[FVGridGeometry::cellCenterIdx()].resize(fvGridGeometry->numCellCenterDofs());
    writeFrom[FVGridGeometry::faceIdx()].resize(fvGridGeometry->numFaceDofs());

    SolutionVector readTo = writeFrom;

    // the problem (initial and boundary conditions)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    assignValues(writeFrom[FVGridGeometry::cellCenterIdx()], values);
    assignValues(writeFrom[FVGridGeometry::faceIdx()], std::array<typename GET_PROP_TYPE(TypeTag, Scalar), 1>{1.0});

    problem->updateStaticWallProperties();
    problem->updateDynamicWallProperties(writeFrom);

    // the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(writeFrom);

    // initialize the vtk output module
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, writeFrom, fileName);
    VtkOutputFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0);

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using IOFields = typename GET_PROP_TYPE(TypeTag, IOFields);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    // cc dofs
    loadSolution(readTo[FVGridGeometry::cellCenterIdx()], fileName + "-00000.vtu",
                 createCellCenterPVNameFunction<IOFields, CellCenterPrimaryVariables, ModelTraits, FluidSystem>(),
                 *fvGridGeometry);

    if (verbose)
    {
        std::cout << "reference cc " << std::endl;
        for (const auto& block : writeFrom[FVGridGeometry::cellCenterIdx()])
        std::cout << block << std::endl;

        std::cout << "result cc " << std::endl;
        for (const auto& block : readTo[FVGridGeometry::cellCenterIdx()])
        std::cout << block << std::endl;
    }

    for (int i = 0; i < readTo[FVGridGeometry::cellCenterIdx()].size(); ++i)
    {
        if (Dune::FloatCmp::ne(readTo[FVGridGeometry::cellCenterIdx()][i], writeFrom[FVGridGeometry::cellCenterIdx()][i]))
            DUNE_THROW(Dune::IOError, "Values don't match: new " << readTo[FVGridGeometry::cellCenterIdx()][i] << ", old " << writeFrom[FVGridGeometry::cellCenterIdx()][i]);
    }

    // face dofs
    loadSolution(readTo[FVGridGeometry::faceIdx()], fileName + "-face-00000.vtp",
                 createFacePVNameFunction<IOFields, FacePrimaryVariables, ModelTraits, FluidSystem>(),
                 *fvGridGeometry);

     if (verbose)
     {
         std::cout << "reference face " << std::endl;
         for (const auto& block : writeFrom[FVGridGeometry::faceIdx()])
         std::cout << block << std::endl;

         std::cout << "result face " << std::endl;
         for (const auto& block : readTo[FVGridGeometry::faceIdx()])
         std::cout << block << std::endl;
     }

     for (int i = 0; i < readTo[FVGridGeometry::faceIdx()].size(); ++i)
     {
         if (Dune::FloatCmp::ne(readTo[FVGridGeometry::faceIdx()][i], writeFrom[FVGridGeometry::faceIdx()][i]))
             DUNE_THROW(Dune::IOError, "Values don't match: new " << readTo[FVGridGeometry::faceIdx()][i] << ", old " << writeFrom[FVGridGeometry::faceIdx()][i]);
     }

     // clean up the folder
     if(deleteFiles)
     {
         if(system("exec rm -r *pvd *vtu *vtp"))
             DUNE_THROW(Dune::IOError, "Deleting files failed");
     }
}


int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize ther parameters
    auto parameters = [](auto& params)
    {
        params["Problem.Name"] = "test";
        params["Vtk.WriteFaceData"] = "true";
    };

    Parameters::init(parameters);

    using CommonTypeTag = TTAG(StaggeredPVNamesTestTypeTag);
    using Grid = typename GET_PROP_TYPE(CommonTypeTag, Grid);
    using FVGridGeometry = typename GET_PROP_TYPE(CommonTypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(CommonTypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, Grid::dimension>;

    const GlobalPosition lowerLeft(0.0);
    const GlobalPosition upperRight(5.0);
    std::array<unsigned int, Grid::dimension> cells;
    std::fill(cells.begin(), cells.end(), 5);

    const auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cells);
    const auto gridView = grid->leafGridView();
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(gridView);
    fvGridGeometry->update();

    using FluidSystem = typename GET_PROP_TYPE(CommonTypeTag, FluidSystem);
    FluidSystem::init();

    testWriteAndReadVtk<TTAG(NavierStokesPVNameTypeTag)>(fvGridGeometry, std::array<Scalar, 1>{1e5}, "navierstokes");
    testWriteAndReadVtk<TTAG(NavierStokesNIPVNameTypeTag)>(fvGridGeometry, std::array<Scalar, 2>{1e5, 300.0}, "navierstokesni");
    testWriteAndReadVtk<TTAG(NavierStokesNCPVNameTypeTag)>(fvGridGeometry, std::array<Scalar, 2>{1e5, 1e-3}, "navierstokesnc");
    testWriteAndReadVtk<TTAG(NavierStokesNCNIPVNameTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1e-3, 300.0}, "navierstokesncni");

    testWriteAndReadVtk<TTAG(ZeroEqNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 1>{1e5}, "zeroeq");
    testWriteAndReadVtk<TTAG(ZeroEqNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 2>{1e5, 300.0}, "zeroeqni");
    testWriteAndReadVtk<TTAG(ZeroEqNCNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 2>{1e5, 1e-3}, "zeroeqnc");
    testWriteAndReadVtk<TTAG(ZeroEqNCNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1e-3, 300.0}, "zeroeqncni");

    testWriteAndReadVtk<TTAG(OneEqNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 2>{1e5, 1.0}, "oneeq");
    testWriteAndReadVtk<TTAG(OneEqNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1.0, 300.0}, "oneeqni");
    testWriteAndReadVtk<TTAG(OneEqNCNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1e-3, 1.0}, "oneeqnc");
    testWriteAndReadVtk<TTAG(OneEqNCNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1e-3, 1.0, 300.0}, "oneeqncni");

    testWriteAndReadVtk<TTAG(KEpsilonNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1.1, 1.2}, "kepsilon");
    testWriteAndReadVtk<TTAG(KEpsilonNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1.1, 1.2, 300.0}, "kepsilonni");
    testWriteAndReadVtk<TTAG(KEpsilonNCNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1e-3, 1.1, 1.2}, "kepsilonnc");
    testWriteAndReadVtk<TTAG(KEpsilonNCNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 5>{1e5, 1e-3, 1.1, 1.2, 300.0}, "kepsilonncni");

    testWriteAndReadVtk<TTAG(LowReKEpsilonNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1.1, 1.2}, "lowrekepsilon");
    testWriteAndReadVtk<TTAG(LowReKEpsilonNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1.1, 1.2, 300.0}, "lowrekepsilonni");
    testWriteAndReadVtk<TTAG(LowReKEpsilonNCNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1e-3, 1.1, 1.2}, "lowrekepsilonnc");
    testWriteAndReadVtk<TTAG(LowReKEpsilonNCNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 5>{1e5, 1e-3, 1.1, 1.2, 300.0}, "lowrekepsilonncni");

    testWriteAndReadVtk<TTAG(KOmegaNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 3>{1e5, 1.1, 1.2}, "komega");
    testWriteAndReadVtk<TTAG(KOmegaNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1.1, 1.2, 300.0}, "komegani");
    testWriteAndReadVtk<TTAG(KOmegaNCNameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 4>{1e5, 1e-3, 1.1, 1.2}, "komeganc");
    testWriteAndReadVtk<TTAG(KOmegaNCNINameTestTypeTag)>(fvGridGeometry, std::array<Scalar, 5>{1e5, 1e-3, 1.1, 1.2, 300.0}, "komegancni");

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 2;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 3;
}
