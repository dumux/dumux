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
 * \brief Test for the staggered grid multi-component RANS model
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

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
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/loadsolution.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>

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
SET_TYPE_PROP(StaggeredPVNamesTestTypeTag, Grid, Dune::YaspGrid<DIM>);

SET_PROP(StaggeredPVNamesTestTypeTag, Problem)
{
    template<class TTag>
    class MockProblem : public NavierStokesProblem<TTag>
    {
        using ParentType = NavierStokesProblem<TTag>;
        using BoundaryTypes = typename GET_PROP_TYPE(TTag, BoundaryTypes);
        using Scalar = typename GET_PROP_TYPE(TTag, Scalar);
    public:
        using ParentType::ParentType;

        template<class GlobalPosition>
        BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
        {
             BoundaryTypes values;
             return values;
        }

        Scalar temperature() const
        {
            return 300;
        }

        void updateStaticWallProperties() {}

        template<class T>
        void updateDynamicWallProperties(const T&) {}
    };

    template<class TTag>
    class TurbulentMockProblem : public KEpsilonProblem<TTag>
    {
        using ParentType = KEpsilonProblem<TTag>;
        using BoundaryTypes = typename GET_PROP_TYPE(TTag, BoundaryTypes);
        using Scalar = typename GET_PROP_TYPE(TTag, Scalar);
    public:
        using ParentType::ParentType;

        template<class GlobalPosition>
        BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
        {
             BoundaryTypes values;
             return values;
        }

        Scalar temperature() const
        {
            return 300;
        }

        template<class Scvf>
        bool isOnWall(const Scvf&) const
        { return true; }
    };

    using Traits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using type = std::conditional_t<Traits::usesTurbulenceModel(), TurbulentMockProblem<TypeTag>, MockProblem<TypeTag>>;
};

}
}

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\nPlease use the provided input files.\n";
        std::cout << errorMessageOut
                  << "\n";
    }
}

template<class ModelTraits, class FluidSystem, class ExpectedResult>
void testNames(const ExpectedResult& expectedResult, bool verbose = true)
{
    // test all primary variable names
    std::cout << "Testing " << Dune::className<ModelTraits>() << std::endl;
    std::cout << "testing all primary variable names " << std::endl;
    for (int pvIdx = 0; pvIdx < ModelTraits::numEq(); ++pvIdx)
    {
        const std::string name = ModelTraits::template primaryVariableName<FluidSystem>(pvIdx);
        if (verbose)
            std::cout << "pvIdx " << pvIdx  << ": " << name << std::endl;
        if (name != expectedResult[pvIdx])
            DUNE_THROW(Dune::IOError, "Wrong primary variable name. Expected " << expectedResult[pvIdx] << ", got " << name);
    }

    // the the names of the primary variables living on the cell centers
    std::cout << "\ntesting cell center primary variable names " << std::endl;
    auto ccPVNames = Dumux::createCellCenterPVNameFunction<ModelTraits, FluidSystem>();
    const auto numCCEq = ModelTraits::numEq() - ModelTraits::dim();
    for (int pvIdx = 0; pvIdx < numCCEq; ++pvIdx)
    {
        const std::string name = ccPVNames(pvIdx, /*state*/0);
        if (verbose)
            std::cout << "pvIdx " << pvIdx  << ": " << name << std::endl;
        if (name != expectedResult[pvIdx + ModelTraits::dim()])
            DUNE_THROW(Dune::IOError, "Wrong primary variable name. Expected " << expectedResult[pvIdx + ModelTraits::dim()] << ", got " << name);
    }

    // the names of the primary variables living on the faces
    std::cout << "\ntesting face primary variable names " << std::endl;
    auto facePVNames = Dumux::createFacePVNameFunction<ModelTraits, FluidSystem>();
    const auto numCCFace = 1;
    for (int pvIdx = 0; pvIdx < numCCFace; ++pvIdx)
    {
        const std::string name = facePVNames(pvIdx, /*state*/0);
        if (verbose)
            std::cout << "pvIdx " << pvIdx  << ": " << name << std::endl;
        if (name != "v")
            DUNE_THROW(Dune::IOError, "Wrong primary variable name. Expected v, got " << name);
    }

    std::cout << std::endl;
}

template<class TypeTag>
void testModel(const std::vector<std::string>&& ccPriVarNames)
{
    auto velocities = [](const int dim)
    {
        if (dim == 1)
            return std::vector<std::string>{"v_x"};
        else if (dim == 2)
            return std::vector<std::string>{"v_x", "v_y"};
        else
            return std::vector<std::string>{"v_x", "v_y", "v_z"};
    };

    // test the methods returning the primary variable names
    using Traits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    auto expectedResult = velocities(Traits::dim());
    expectedResult.insert(expectedResult.end(), ccPriVarNames.begin(), ccPriVarNames.end());

    testNames<Traits, FluidSystem>(expectedResult);
}

template<class SolutionVector, class Values>
void assignValues(SolutionVector& sol, Values values)
{
    for (auto& entry : sol)
    {
        for (int pvIdx = 0; pvIdx < decltype(values)::dimension; ++pvIdx)
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

template<class TypeTag, class FVGridGeometry, class Values>
void testWriteAndReadVtk(std::shared_ptr<FVGridGeometry> fvGridGeometry, const Values& values, const std::string& fileName, bool verbose = false)
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

    problem->updateStaticWallProperties();
    problem->updateDynamicWallProperties(writeFrom);

    // the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(writeFrom);

    const std::string fileNameCell = fileName + "_cell";

    // initialize the vtk output module
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, writeFrom, fileNameCell);
    VtkOutputFields::init(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0);

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    loadSolution(readTo[FVGridGeometry::cellCenterIdx()], fileNameCell + "-00000.vtu",
                 createCellCenterPVNameFunction<ModelTraits, FluidSystem>(),
                 *fvGridGeometry);

    if (verbose)
    {
        std::cout << "old " << std::endl;
        for (const auto& block : writeFrom[FVGridGeometry::cellCenterIdx()])
        std::cout << block << std::endl;

        std::cout << "new " << std::endl;
        for (const auto& block : readTo[FVGridGeometry::cellCenterIdx()])
        std::cout << block << std::endl;
    }

    for (int i = 0; i < readTo[FVGridGeometry::cellCenterIdx()].size(); ++i)
    {
        if (Dune::FloatCmp::ne(readTo[FVGridGeometry::cellCenterIdx()][i], writeFrom[FVGridGeometry::cellCenterIdx()][i]))
            DUNE_THROW(Dune::IOError, "Values don't match: new " << readTo[FVGridGeometry::cellCenterIdx()][i] << ", old " << writeFrom[FVGridGeometry::cellCenterIdx()][i]);
    }

    // TODO face dofs

}


int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

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

    {
        using TypeTag = TTAG(NavierStokesPVNameTypeTag);
        testModel<TypeTag>({"p"});
        Dune::FieldVector<Scalar, 1> values = {1e5};
        testWriteAndReadVtk<TypeTag>(fvGridGeometry, values, "navierstokes");
    }

    {
        using TypeTag = TTAG(NavierStokesNIPVNameTypeTag);
        testModel<TypeTag>({"p", "T"});
    }

    {
        using TypeTag = TTAG(NavierStokesNCPVNameTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas"});
    }

    {
        using TypeTag = TTAG(NavierStokesNCNIPVNameTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "T"});
    }

    // KEpsilon
    {
        using TypeTag = TTAG(KEpsilonNameTestTypeTag);
        testModel<TypeTag>({"p", "k", "epsilon"});
    }

    {
        using TypeTag = TTAG(KEpsilonNINameTestTypeTag);
        testModel<TypeTag>({"p", "k", "epsilon", "T"});
    }

    {
        using TypeTag = TTAG(KEpsilonNCNameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "k", "epsilon"});
    }

    {
        using TypeTag = TTAG(KEpsilonNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "k", "epsilon", "T"});
        Dune::FieldVector<Scalar, 5> values = {1e5, 1e-3, 1.0, 1.0, 300};
        testWriteAndReadVtk<TypeTag>(fvGridGeometry, values, "navierstokes");
    }

    // Low-Re-KEpsilon
    {
        using TypeTag = TTAG(LowReKEpsilonNameTestTypeTag);
        testModel<TypeTag>({"p", "k", "epsilon"});
    }

    {
        using TypeTag = TTAG(LowReKEpsilonNINameTestTypeTag);
        testModel<TypeTag>({"p", "k", "epsilon", "T"});
    }

    {
        using TypeTag = TTAG(LowReKEpsilonNCNameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "k", "epsilon"});
    }

    {
        using TypeTag = TTAG(LowReKEpsilonNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "k", "epsilon", "T"});
    }

    // KOmega
    {
        using TypeTag = TTAG(KOmegaNameTestTypeTag);
        testModel<TypeTag>({"p", "k", "omega"});
    }

    {
        using TypeTag = TTAG(KOmegaNINameTestTypeTag);
        testModel<TypeTag>({"p", "k", "omega", "T"});
    }

    {
        using TypeTag = TTAG(KOmegaNCNameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "k", "omega"});
    }

    {
        using TypeTag = TTAG(KOmegaNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "k", "omega", "T"});
    }

    // One-Eq
    {
        using TypeTag = TTAG(OneEqNameTestTypeTag);
        testModel<TypeTag>({"p", "nu_tilde"});
    }

    {
        using TypeTag = TTAG(OneEqNINameTestTypeTag);
        testModel<TypeTag>({"p", "nu_tilde", "T"});
    }

    {
        using TypeTag = TTAG(OneEqNCNameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "nu_tilde"});
    }

    {
        using TypeTag = TTAG(OneEqNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "nu_tilde", "T"});
    }

    // Zero-Eq
    {
        using TypeTag = TTAG(ZeroEqNameTestTypeTag);
        testModel<TypeTag>({"p"});
    }

    {
        using TypeTag = TTAG(ZeroEqNINameTestTypeTag);
        testModel<TypeTag>({"p","T"});
    }

    {
        using TypeTag = TTAG(ZeroEqNCNameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas"});
    }

    {
        using TypeTag = TTAG(ZeroEqNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "X^H2O_gas", "T"});
    }

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
