import argparse
import os, sys
import subprocess
from fuzzycomparevtu import compare_vtk

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--command', nargs=1, help='The executable and optional arguments as a single string', required=True)
parser.add_argument('-s', '--script', nargs=1, help="The comparison script. [fuzzy, exact, <path_to_script>] where the script takes two vtu files as arguments.")
parser.add_argument('-f', '--files', nargs='+', help="Pairs of reference and vtu file names. Usage: '[-f ref1 vtu1 [[ref2] [vtu2] ...]]'")
parser.add_argument('-r', '--relative', type=float, default=1e-2, help='maximum relative error (default=1e-2) when using fuzzy comparison')
parser.add_argument('-a', '--absolute', type=float, default=1e-9, help='maximum absolute error (default=1e-9) when using fuzzy comparison')
args = vars(parser.parse_args())

# check parameters
if args['script']:
    if len(args['files'])%2 != 0 or not args['files']:
        sys.stderr.write("The files have to be pairs of vtu and reference files. Usage '-f [ref1] [vtu1] [[ref2] [vtu2] ...]'")
        parser.print_help()
        sys.exit(1)
    for i in range(0, len(args['files'])//2):
        # delete the vtu files to compare
        ref_dir = os.path.dirname(os.path.abspath(__file__)).rstrip("bin") + "test/references"
        if os.path.dirname(args['files'][(i*2)+1]) == ref_dir:
            sys.stderr.write("Tried to delete a reference solution. Specify reference file first, then the VTU file. Usage: '[-f ref1 vtu1 [[ref2] [vtu2] ...]]'")
            sys.exit(1)
        subprocess.call(['rm', '-fv', args['files'][(i*2)+1]])

# run the test
res = 1
try:
    res = subprocess.call(args['command'][0].split())
except OSError as e:
    print(args['command'][0].split())
    print("OSError: Command not found. Most likely the executable specified doesn't exist.")
    sys.exit(1)
if res:
    sys.exit(1)

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
            if args['relative'] and args['absolute']:
                result = compare_vtk(args['files'][i*2], args['files'][(i*2)+1], relative=args['relative'], absolute=args['absolute'])
            elif args['relative'] and not args['absolute']:
                result = compare_vtk(args['files'][i*2], args['files'][(i*2)+1], relative=args['relative'])
            elif args['absolute'] and not args['relative']:
                result = compare_vtk(args['files'][i*2], args['files'][(i*2)+1], absolute=args['absolute'])
            else:
                result = compare_vtk(args['files'][i*2], args['files'][(i*2)+1])
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
