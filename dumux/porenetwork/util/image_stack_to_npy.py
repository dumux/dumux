#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Reads a stack of binarized images representing slices of a 3D domain and stores the data
in a binary 3D numpy array.
"""

from PIL import Image
import numpy as np
import os.path
import sys
import argparse


def read_images(filenames, invert=False):
    arr3d = []
    for filename in filenames:
        if not os.path.isfile(filename):
            print(f"\rimage {filename} not found")
            continue
        message = f"reading image {filename}"
        print(f"\r{message}", end="")
        im = Image.open(filename)
        a = np.asarray(im)
        im.close()
        arr3d.append(a)
        print(f"\r{' '*len(message)}", end="")
    arr3d = np.asarray(arr3d) / 255
    if invert:
        arr3d = 1 - arr3d
    print(f"\rRead {arr3d.shape[0]}/{len(filenames)} images.")
    return arr3d


def read_image_stack(format_string, n_start, n_end, invert=False):
    filenames = [format_string.format(i) for i in range(n_start, n_end)]
    return read_images(filenames, invert)


def main():
    parser = argparse.ArgumentParser(
        description="Process a stack of images into a binary numpy array"
    )
    parser.add_argument("outfile", help="filename of output file")
    parser.add_argument(
        "filenames", metavar="files", nargs="*", help="filenames of images to include"
    )
    parser.add_argument(
        "--format",
        dest="format_string",
        nargs="?",
        help="format string with one unnamed field for slice index",
    )
    parser.add_argument(
        "--from", dest="start", type=int, nargs="?", default=0, help="starting index"
    )
    parser.add_argument(
        "--to", dest="end", type=int, nargs="?", default=0, help="ending index (exclusive)"
    )
    parser.add_argument(
        "--invert", dest="invert", action="store_true", default=False, help="invert pore and solid"
    )
    args = parser.parse_args()
    filenames = args.filenames
    if args.format_string:
        filenames += [args.format_string.format(i) for i in range(args.start, args.end)]
    array3d = read_images(filenames, args.invert)
    np.save(args.outfile, array3d)


if __name__ == "__main__":
    main()
