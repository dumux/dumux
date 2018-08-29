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

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
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

namespace Dumux {
namespace Properties {

NEW_TYPE_TAG(StaggeredPVNamesTestTypeTag);

NEW_TYPE_TAG(NavierStokesPVNameTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokes, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(NavierStokesNIPVNameTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNI, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(NavierStokesNCPVNameTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNC, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(NavierStokesNCNIPVNameTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNCNI, StaggeredPVNamesTestTypeTag));

NEW_TYPE_TAG(KEpsilonNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilon, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(KEpsilonNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilonNI, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(KEpsilonNCNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilonNC, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(KEpsilonNCNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilonNCNI, StaggeredPVNamesTestTypeTag));

NEW_TYPE_TAG(LowReKEpsilonNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilon, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(LowReKEpsilonNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilonNI, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(LowReKEpsilonNCNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilonNC, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(LowReKEpsilonNCNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilonNCNI, StaggeredPVNamesTestTypeTag));

NEW_TYPE_TAG(KOmegaNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KOmega, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(KOmegaNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KOmegaNI, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(KOmegaNCNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KOmegaNC, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(KOmegaNCNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KOmegaNCNI, StaggeredPVNamesTestTypeTag));

NEW_TYPE_TAG(OneEqNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEq, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(OneEqNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNI, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(OneEqNCNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNC, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(OneEqNCNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNCNI, StaggeredPVNamesTestTypeTag));

NEW_TYPE_TAG(ZeroEqNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEq, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(ZeroEqNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNI, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(ZeroEqNCNameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNC, StaggeredPVNamesTestTypeTag));
NEW_TYPE_TAG(ZeroEqNCNINameTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNCNI, StaggeredPVNamesTestTypeTag));

// The fluid system
SET_PROP(StaggeredPVNamesTestTypeTag, FluidSystem)
{
  using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
SET_TYPE_PROP(StaggeredPVNamesTestTypeTag, Grid, Dune::YaspGrid<DIM>);
SET_BOOL_PROP(StaggeredPVNamesTestTypeTag, UseMoles, true);

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

    // define the type tag for this problem
    // using TypeTag = TTAG(KEpsilonNCNameTestTypeTag);

    // test the methods returning the primary variable names
    using Traits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    auto expectedResult = velocities(Traits::dim());
    expectedResult.insert(expectedResult.end(), ccPriVarNames.begin(), ccPriVarNames.end());

    testNames<Traits, FluidSystem>(expectedResult);
}


int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    {
        using TypeTag = TTAG(NavierStokesPVNameTypeTag);
        testModel<TypeTag>({"p"});
    }

    {
        using TypeTag = TTAG(NavierStokesNIPVNameTypeTag);
        testModel<TypeTag>({"p", "T"});
    }

    {
        using TypeTag = TTAG(NavierStokesNCPVNameTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas"});
    }

    {
        using TypeTag = TTAG(NavierStokesNCNIPVNameTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas", "T"});
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
        testModel<TypeTag>({"p", "x^H2O_gas", "k", "epsilon"});
    }

    {
        using TypeTag = TTAG(KEpsilonNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas", "k", "epsilon", "T"});
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
        testModel<TypeTag>({"p", "x^H2O_gas", "k", "epsilon"});
    }

    {
        using TypeTag = TTAG(LowReKEpsilonNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas", "k", "epsilon", "T"});
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
        testModel<TypeTag>({"p", "x^H2O_gas", "k", "omega"});
    }

    {
        using TypeTag = TTAG(KOmegaNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas", "k", "omega", "T"});
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
        testModel<TypeTag>({"p", "x^H2O_gas", "nu_tilde"});
    }

    {
        using TypeTag = TTAG(OneEqNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas", "nu_tilde", "T"});
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
        testModel<TypeTag>({"p", "x^H2O_gas"});
    }

    {
        using TypeTag = TTAG(ZeroEqNCNINameTestTypeTag);
        testModel<TypeTag>({"p", "x^H2O_gas", "T"});
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
