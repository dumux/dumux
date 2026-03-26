// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4 sw=4 sts=4
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Utility functions for simple parameter file parsing (JSON-like syntax).
 */
// jsonparametertools.hh - Converts nested JSON to flat Dune params!
#ifndef DUMUX_JSONPARAMETERTOOLS_HH
#define DUMUX_JSONPARAMETERTOOLS_HH

#include <dune/common/parametertree.hh>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>

namespace Dumux {
namespace Utils {

class JsonToTreeConverter {
    std::string data;
    size_t pos = 0;

    void skipWhitespace() {
        while (pos < data.size() && (std::isspace(data[pos]) || (unsigned char)data[pos] == 0xEF || (unsigned char)data[pos] == 0xBB || (unsigned char)data[pos] == 0xBF)) pos++;
    }

    std::string parseString() {
        skipWhitespace();
        if (pos < data.size() && data[pos] == '\"') pos++;
        std::string s;
        while (pos < data.size() && data[pos] != '\"') {
            s += data[pos++];
        }
        if (pos < data.size() && data[pos] == '\"') pos++;
        return s;
    }

    void parseValue(const std::string& prefix, Dune::ParameterTree& pt) {
        skipWhitespace();
        if (pos >= data.size()) return;

        if (data[pos] == '{') {
            pos++;
            while (pos < data.size() && data[pos] != '}') {
                std::string key = parseString();
                skipWhitespace();
                if (pos < data.size() && data[pos] == ':') pos++;

                std::string nextPrefix = prefix.empty() ? key : prefix + "." + key;
                parseValue(nextPrefix, pt);

                skipWhitespace();
                if (pos < data.size() && data[pos] == ',') pos++;
                skipWhitespace();
            }
            if (pos < data.size() && data[pos] == '}') pos++;
        } else {
            std::string val;
            if (data[pos] == '\"') {
                val = parseString();
            } else {
                while (pos < data.size() && data[pos] != ',' && data[pos] != '}' && data[pos] != ']') {
                    if (!std::isspace(data[pos])) val += data[pos];
                    pos++;
                }
            }
            if (!prefix.empty()) pt[prefix] = val;
        }
    }

public:
    void convert(const std::string& filename, Dune::ParameterTree& pt) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "!!! FILE ERROR: Could not open " << filename << ". Check if it exists in " << std::string(getenv("PWD")) << std::endl;
            return;
        }

        data.assign((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        file.close();

        if (data.empty()) {
            std::cerr << "!!! FILE ERROR: File " << filename << " is empty!" << std::endl;
            return;
        }

        pos = 0;
        parseValue("", pt);
    }
};

inline void parseNestedJson(const std::string& filename, Dune::ParameterTree& pt) {
    JsonToTreeConverter converter;
    converter.convert(filename, pt);
}

} // Utils
} // Dumux
#endif