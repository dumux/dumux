// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief This file tests the metadata functionality
 */
#include <config.h>

#include <string>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/metadata.hh>

struct MyDataClass
{
    double val;
    std::string name;
};

int main()
{
    using namespace Dumux::MetaData;

    Collector collector;

    // Test brackets operator and Json tree
    collector["MyData"]["values"] = {1,2};
    collector["MyData"]["myBool"] = true;

    if(!collector["MyData"]["myBool"])
        DUNE_THROW(Dune::Exception, "Wrong bool in JsonTree");

    auto& vals = collector["MyData"]["values"];
    vals.push_back(3);

    for(int i=0; i<vals.size(); i++)
        if (vals[i] != i+1)
            DUNE_THROW(Dune::Exception, "Wrong vector in JsonTree");

    // Test className function
    MyDataClass data({3.1,"mydata"});
    std::string className("MyDataClass");
    collector[className]["className"] = Collector::className(data, true);
    collector[className]["value"] = data.val;
    collector[className]["name"] = data.name;

    if (!(collector[className]["className"].get<std::string>().find(className) != std::string::npos))
        DUNE_THROW(Dune::Exception, "Class name is wrong");

    if (!Dune::FloatCmp::eq(3.1, collector[className]["value"].get<double>() , 1e-6))
        DUNE_THROW(Dune::Exception, "Wrong double value in JsonTree");

    std::string s(collector[className]["name"].get<std::string>());
    if(!(data.name.compare(s)==0))
        DUNE_THROW(Dune::Exception, "Wrong string in JsonTree");

    // Test writing and reading Json files
    writeJsonFile(collector, "output");
    Collector collector2;
    readJsonFile(collector2, "output");

    if(!(collector.getTree() == collector2.getTree()))
        DUNE_THROW(Dune::Exception, "Json trees differ");

    collector2[className]["value"] = 3.2;
    if((collector.getTree() == collector2.getTree()))
        DUNE_THROW(Dune::Exception, "Json trees do not differ");

    return 0;
}
