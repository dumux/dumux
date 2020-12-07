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
 * \brief Test writing and reading sequence container to and from file
 */
#include <config.h>
#include <vector>
#include <algorithm>
#include <iostream>

#include <dumux/io/rasterimagereader.hh>
#include <dumux/io/container.hh>

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    using namespace Dumux;

    //////////////////////////////////////////////////
    // Test the binary (black and white) image reader
    //////////////////////////////////////////////////

    // read an ASCII image without applying Dune's cell ordering (origin at upper left)
    auto blackAndWhiteImageASCII = NetPBMReader::readPBM("blackwhite_j.pbm", false);

    // create a 2D array to show the image in the terminal
    std::vector<std::vector<bool>> printableBlackAndWhiteImage;
    printableBlackAndWhiteImage.resize(blackAndWhiteImageASCII.header().nRows, std::vector<bool>(blackAndWhiteImageASCII.header().nCols));
    NetPBMReader::fillImage(printableBlackAndWhiteImage, blackAndWhiteImageASCII);
    // print the image
    for (const auto& row: printableBlackAndWhiteImage)
    {
        for (const auto& col : row)
            std::cout << col << " ";
        std::cout << std::endl;
    }

    // helper function to determine if the data of two vector-like classes are equal
    auto isEqual = [](const auto& v1, const auto& v2)
    {
        return std::equal(v1.begin(), v1.end(), v2.begin());
    };

    // test the flattenImageToVector helper function
    if (!isEqual(blackAndWhiteImageASCII, NetPBMReader::flattenImageToVector(printableBlackAndWhiteImage)))
    {
        std::cout << "flattenImageToVector failed" << std::endl;
        return 1;
    }

    std::cout << std::endl;

    // change the ordering according to Dune's convention (origin at lower left)
    NetPBMReader::applyDuneGridOrdering(blackAndWhiteImageASCII);
    // compare the result to a reference stored in a text file
    if (!isEqual(blackAndWhiteImageASCII, readFileToContainer<std::vector<bool>>("blackwhite_j.txt")))
    {
        std::cout << "Reading black/white (ASCII) failed" << std::endl;
        return 1;
    }

    std::cout << std::endl;

    // read a binary image and directly apply Dune's ordering
    const auto blackAndWhiteImageBinary = NetPBMReader::readPBM("blackwhite_binary_j.pbm");
    if (!isEqual(blackAndWhiteImageBinary, blackAndWhiteImageASCII))
    {
        std::cout << "Reading black/white (binary) failed" << std::endl;
        return 1;
    }

    std::cout << std::endl;

    // test file where the dimensions are given in the same line as the magic number and comments are present
    const std::vector<bool> reference{1,0,0,1};

    if (!isEqual(reference, NetPBMReader::readPBM("blackwhite_dim_firstline.pbm", false)))
    {
        std::cout << "Reading black/white with dimension in first line failed" << std::endl;
        return 1;
    }
    if (!isEqual(reference, NetPBMReader::readPBM("blackwhite_dim_firstline.pbm", false)))
    {
        std::cout << "Reading black/white (binary) with dimension in first line failed" << std::endl;
        return 1;
    }

    // test error message for poorly formatted file
    try
    {
        NetPBMReader::readPBM("blackwhite_fail.pbm", false);
    }
    catch(const Dune::IOError& e)
    {
        const auto tokens = tokenize(e.what(), "]:");
        if (tokens.back() != " Expecting only dimensions (2 numbers) in line 3")
        {
            std::cout << e.what() << std::endl;
            return 1;
        }
    }

    //////////////////////////////////////////////////
    // Test the gray scale image reader
    //////////////////////////////////////////////////

    // read an ASCII image without applying Dune's cell ordering (origin at upper left)
    auto grayScaleImageASCII = NetPBMReader::template readPGM<std::size_t>("grayscale_j.pgm", false);
    std::vector<std::vector<std::uint8_t>> printableGrayScaleImage;
    printableGrayScaleImage.resize(grayScaleImageASCII.header().nRows, std::vector<std::uint8_t>(grayScaleImageASCII.header().nCols));
    NetPBMReader::fillImage(printableGrayScaleImage, grayScaleImageASCII);
    // print the image
    for (const auto& row: printableGrayScaleImage)
    {
        for (const auto& col : row)
            std::cout << std::setw(3) << +col << " ";
        std::cout << std::endl;
    }

    // change the ordering according to Dune's convention (origin at lower left)
    NetPBMReader::applyDuneGridOrdering(grayScaleImageASCII);
    // compare the result to a reference stored in a text file
    if (!isEqual(grayScaleImageASCII, readFileToContainer<std::vector<std::size_t>>("grayscale_j.txt")))
    {
        std::cout << "Reading black/white (ASCII) failed" << std::endl;
        return 1;
    }

    std::cout << std::endl;

    // read a binary image and directly apply Dune's ordering
    const auto grayScaleImageBinary = NetPBMReader::readPGM<std::size_t>("grayscale_binary_j.pgm");
    if (!isEqual(grayScaleImageBinary, grayScaleImageASCII))
    {
        std::cout << "Reading black/white (binary) failed" << std::endl;
        return 1;
    }

    // test error message for poorly formatted file
    try
    {
        NetPBMReader::readPBM("grayscale_fail_binary.pgm", false);
    }
    catch(const Dune::IOError& e)
    {
        const auto tokens = tokenize(e.what(), "]:");
        if (tokens.back() != " Expecting only intensity (one number) in line 4")
        {
            std::cout << e.what() << std::endl;
            return 1;
        }
    }

    return 0;
}
