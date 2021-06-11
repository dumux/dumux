#!/usr/bin/env python3

"""
A simple script for conveniently running the SNOW algorithm of PoreSpy
(watershed image segmentation) in order to extract a (dual) pore network
from a 2D or 3D raster image. Accepts .raw and .nrrd images.
Requires the PoreSpy (https://porespy.org/) and OpenPNM package (https://openpnm.org/).
"""

from argparse import ArgumentParser
import numpy as np
import openpnm as op
import porespy as ps

parser = ArgumentParser(description='extract network from image using PoreSpy')
parser.add_argument('file', type=str, help='The image file')
parser.add_argument('-r', '--resolution', type=float, help='The pixel/voxel size', required=True)
parser.add_argument('-s', '--size', nargs='+', type=int, help='The number of pixels/voxels in each direction', default=[])
parser.add_argument('-n', '--outputname', type=str, help='The output file name', default='')
parser.add_argument('-m', '--marchingCubesArea', action='store_true', help="Use marching cubes alg. to calculate throat area")
parser.add_argument('-d', '--dualSnow', action='store_true', help="Extract dual network")
parser.add_argument('-i', '--invert', action='store_true', help="Invert the image data")
args = vars(parser.parse_args())

if args['file'].endswith('raw'):
    if not args['outputname']:
        args['outputname'] = args['file'].replace('.raw', '')

    if len(args['size']) < 2 or len(args['size']) > 3:
        raise IOError('argument --voxels requires numVoxels per dimension')

    im = np.fromfile(args['file'], dtype='uint8', sep="").reshape(args['size'])
    im = np.array(im, dtype=bool)

    if len(args['size']) == 2:
        dim = 2
        print("Read 2D raw image file", args['file'], "with", args['size'], "pixels")
    else:
        dim = 3
        print("Read 3D raw image file", args['file'], "with", args['size'], "voxels")
        im = np.swapaxes(im, 0, 2)

elif args['file'].endswith('nrrd'):
    import nrrd
    im, header = nrrd.read(args['file'])
    im = ~np.array(im, dtype=bool)
    args['size'] = header['sizes']
    print("Reading nrrd image file", args['file'], "with", args['size'], "voxels")
    if not args['outputname']:
        args['outputname'] = args['file'].replace('.nrrd', '')

else:
    raise IOError('Invalid file format')

if args['invert']:
    print('Inverting data')
    im = ~im

if dim == 2:
    ps.io.to_vtk(np.array(im, dtype=int)[:, :, np.newaxis],
                 args['outputname'], voxel_size=args['resolution'])

print("Porosity is", ps.metrics.porosity(im))

if args['marchingCubesArea']:
    print('Using Marching Cubes Alg. to calculate throat areas')

# application of SNOW algorithm
if not args['dualSnow']:
    snowOutput = ps.networks.snow2(im,
                                   voxel_size=args['resolution'],
                                   accuracy=('high' if args['marchingCubesArea'] else 'standard'))
else:
    snowOutput = ps.networks.snow2(im+1,
                                   voxel_size=args['resolution'],
                                   accuracy=('high' if args['marchingCubesArea'] else 'standard'),
                                   phase_alias={1: 'solid', 2: 'void'})

# use PoreSpy to sanitize some parameters (might become obsolete as some point)
pn, geo = op.io.PoreSpy.import_data(snowOutput.network)

# trimming pore network to avoid singularity
print('Number of pores before trimming: ', pn.Np)
h = pn.check_network_health()
op.topotools.trim(network=pn, pores=h['trim_pores'])
print('Number of pores after trimming: ', pn.Np)

filename = args['outputname'] + ('_dual' if args['dualSnow'] else '')

# export to vtk
op.io.VTK.save(network=pn, filename=filename)
print("Writing", filename + ".vtp")

# export network for OpenPNM
pn.project.save_project(filename)
print("Writing", filename + ".pnm")
