""" A module for fuzzy comparing data files.

This module provides methods to compare two data files.
Applicable for all style formats like e.g. csv files.
Fuzzy compares numbers by using absolute and/or relative difference comparison.

"""
import argparse
import csv
import json
import sys
from fuzzycomparevtu import is_fuzzy_equal_text

def compare_data(dataFile1, dataFile2, delimiter, absolute=1.5e-7, relative=1e-2, zeroValueThreshold={}, verbose=True):
    """ take two data files and compare them. Returns an exit key as returnvalue.

    Arguments:
    ----------
    dataFile1, dataFile2 : string
        The filenames of the data files to compare
    delimiter: string
        The delimiter for the columns

    Keyword Arguments:
    ------------------
    absolute : float
        The epsilon used for comparing numbers with an absolute criterion
    relative: float
        The epsilon used for comparing numbers with an relative criterion
    zeroValueThreshold: dict
        A dictionary of parameter value pairs that set the threshold under
        which a number is treated as zero for a certain parameter. Use this parameter if
        you have to avoid comparisons of very small numbers for a certain parameter.
    verbose : bool
        If the script should produce informative output. Enabled by default as the details
        give the tester a lot more information on why tests fail.
    """

    if verbose:
        print("Comparing {} and {}".format(dataFile1, dataFile2))
        print("... with a maximum relative error of {} and a maximum absolute error of {}*max_abs_parameter_value.".format(relative, absolute))

    # construct element tree from data files
    data1 = list(csv.reader(open(dataFile1, 'r'), delimiter=delimiter))
    data2 = list(csv.reader(open(dataFile2, 'r'), delimiter=delimiter))

    if (len(data1) != len(data2)):
        print("Length of data1 and data2 not equal: ref=", len(data1), ",new=", len(data2), ". Aborting!")
        exit (3)

    is_equal = True
    for i in range(0,len(data1[0])):
        a = data1[0][i]
        b = data2[0][i]
        for j in range(1,len(data1)):
            a += " {0}".format(data1[j][i])
            b += " {0}".format(data2[j][i])

        if not is_fuzzy_equal_text(a, b, "row {0}".format(i), len(data1), absolute, relative, zeroValueThreshold, verbose):
            if verbose:
                is_equal = False
            else:
                return False

    if is_equal:
        return 0
    else:
        return 1


# main program if called as script return appropriate error codes
if __name__ == "__main__":
    # handle arguments and print help message
    parser = argparse.ArgumentParser(description='Fuzzy compare of two data files (e.g csv). \
        The files are accepted if for every value the difference is below the absolute error \
        or below the relative error or below both.')
    parser.add_argument('data_file_1', type=str, help='first file to compare')
    parser.add_argument('data_file_2', type=str, help='second file to compare')
    parser.add_argument('delimiter', type=str, help='second file to compare')
    parser.add_argument('-r', '--relative', type=float, default=1e-2, help='maximum relative error (default=1e-2)')
    parser.add_argument('-a', '--absolute', type=float, default=1.5e-7, help='maximum absolute error (default=1.5e-7)')
    parser.add_argument('-v', '--verbose', type=bool, default=True, help='verbosity of the script')
    parser.add_argument('-z', '--zeroThreshold', type=json.loads, default='{}', help='Thresholds for treating numbers as zero for a parameter as a python dict e.g. {"vel":1e-7,"delP":1.0}')
    args = vars(parser.parse_args())

    sys.exit(compare_data(args["data_file_1"], args["data_file_2"], args["delimiter"], args["absolute"], args["relative"], args["zeroThreshold"], args["verbose"]))
