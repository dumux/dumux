// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test writing and reading sequence container to and from file
 */
#include <config.h>
#include <vector>
#include <algorithm>
#include <iostream>

#include <dumux/io/rasterimagereader.hh>
#include <dumux/io/rasterimagewriter.hh>
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
    auto blackAndWhiteImageASCIIOutput = blackAndWhiteImageASCII;
    for (int i = 0; i < blackAndWhiteImageASCII.size(); i++)
    {
        if (i < 5)
            blackAndWhiteImageASCIIOutput[i] = 1;
    }
    NetPBMWriter::write("blackwhite_j_new.pbm", blackAndWhiteImageASCIIOutput, false);

    // create a 2D array to show the image in the terminal
    std::vector<std::vector<bool>> printableBlackAndWhiteImage;
    printableBlackAndWhiteImage.resize(blackAndWhiteImageASCII.header().nRows, std::vector<bool>(blackAndWhiteImageASCII.header().nCols));
    NetPBMReader::fillImage(printableBlackAndWhiteImage, blackAndWhiteImageASCII);
    // print the image
    for (const auto& row: printableBlackAndWhiteImage)
    {
        for (bool col : row)
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
    auto grayScaleImageASCIIOutput = grayScaleImageASCII;
    for (int i = 0; i < grayScaleImageASCII.size(); i++)
    {
        if (i < 5)
            grayScaleImageASCIIOutput[i] = 125;
    }
    NetPBMWriter::write("grayscale_j_new.pgm", grayScaleImageASCIIOutput, false);

    std::vector<std::vector<std::uint8_t>> printableGrayScaleImage;
    printableGrayScaleImage.resize(grayScaleImageASCII.header().nRows, std::vector<std::uint8_t>(grayScaleImageASCII.header().nCols));
    NetPBMReader::fillImage(printableGrayScaleImage, grayScaleImageASCII);
    // print the image
    for (const auto& row: printableGrayScaleImage)
    {
        for (bool col : row)
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
