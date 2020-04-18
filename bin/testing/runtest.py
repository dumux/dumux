#!/usr/bin/env python3
import argparse
import shlex
import os, sys
import subprocess
import json
from fuzzycomparevtu import compare_vtk
from fuzzycomparedata import compare_data

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--command', nargs=1, help='The executable and optional arguments as a single string', required=True)
parser.add_argument('-s', '--script', nargs=1, help="The comparison script. [fuzzy, fuzzyData, exact, <path_to_script>] where the script takes two files as arguments.")
parser.add_argument('-f', '--files', nargs='+', help="Pairs of file names (first reference, then current). Usage: '[-f ref1 cur1 [[ref2] [cur2] ...]]'")
parser.add_argument('-d', '--delimiter', type=str, default=',', help='Column delimiter for data files')
parser.add_argument('-r', '--relative', type=float, default=1e-2, help='maximum relative error (default=1e-2) when using fuzzy comparison')
parser.add_argument('-a', '--absolute', type=float, default=1.5e-7, help='maximum absolute error (default=1.5e-7) when using fuzzy comparison')
parser.add_argument('-z', '--zeroThreshold', type=json.loads, default='{}', help='Thresholds for treating numbers as zero for a parameter as a python dict e.g. {"vel":1e-7,"delP":1.0}')
args = vars(parser.parse_args())

# check parameters
if args['script']:
    if len(args['files'])%2 != 0 or not args['files']:
        sys.stderr.write("The files have to be pairs of reference and current solution files. Usage '-f [ref1] [cur1] [[ref2] [cur2] ...]'")
        parser.print_help()
        sys.exit(1)
    for i in range(0, len(args['files'])//2):
        # delete the vtu files to compare
        ref_dir = os.path.dirname(os.path.abspath(__file__)).rstrip("bin") + "test/references"
        if os.path.dirname(args['files'][(i*2)+1]) == ref_dir:
            sys.stderr.write("Tried to delete a reference solution. Specify reference file first, then the current solution. Usage: '[-f ref1 cur1 [[ref2] [cur2] ...]]'")
            sys.exit(1)
        subprocess.call(['rm', '-fv', args['files'][(i*2)+1]])

# run the test
res = 1
try:
    res = subprocess.call(shlex.split(args['command'][0]))
except OSError as e:
    print(args['command'][0].split())
    print("OSError: Command not found. Most likely the executable specified doesn't exist.")
    sys.exit(1)
if res:
    sys.exit(res)

# run the comparison
if args['script']:
    # exact comparison?
    if args['script'] == ['exact']:
        return_code = 0
        for i in range(0, len(args['files'])//2):
            print("\nExact comparison...")
            result = subprocess.call(['diff', args['files'][i*2], args['files'][(i*2)+1]])
            if result:
                return_code = 1
        sys.exit(return_code)

    # fuzzy comparison?
    elif args['script'] == ["fuzzy"] or args['script'] == [os.path.dirname(os.path.abspath(__file__)) + "/fuzzycomparevtu.py"]:
        return_code = 0
        for i in range(0, len(args['files'])//2):
            print("\nFuzzy comparison...")
            result = compare_vtk(args['files'][i*2], args['files'][(i*2)+1], relative=args['relative'], absolute=args['absolute'], zeroValueThreshold=args['zeroThreshold'])
            if result:
                return_code = 1
        sys.exit(return_code)

    # fuzzy comparison of data sets?
    elif args['script'] == ["fuzzyData"]:
        return_code = 0
        for i in range(0, len(args['files'])//2):
            print("\nFuzzy data comparison...")
            result = compare_data(args['files'][i*2], args['files'][(i*2)+1], args['delimiter'], relative=args['relative'], absolute=args['absolute'], zeroValueThreshold=args['zeroThreshold'])
            if result:
                return_code = 1
        sys.exit(return_code)

    # other script?
    else:
        return_code = 0
        for i in range(0, len(args['files'])//2):
            print("\n{} comparison...".format(args['script']))
            result = subprocess.call(args['script'], args['files'][i*2], args['files'][(i*2)+1])
            if result:
                return_code = 1
        sys.exit(return_code)

# everything is fine
sys.exit(0)
