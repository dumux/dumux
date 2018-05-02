""" A module for fuzzy comparing VTK files.

This module provides methods to compare two VTK files. Applicable
for all VTK style formats like VTK files. Fuzzy compares numbers by
using absolute and/or relative difference comparison.

"""
from __future__ import print_function
from __future__ import absolute_import
import argparse
import xml.etree.ElementTree as ET
from operator import attrgetter, itemgetter
import json
import sys
import math
from six.moves import range
from six.moves import zip


# fuzzy compare VTK tree from VTK strings
def compare_vtk(vtk1, vtk2, absolute=1.5e-7, relative=1e-2, zeroValueThreshold={}, verbose=True):
    """ take two vtk files and compare them. Returns an exit key as returnvalue.

    Arguments:
    ----------
    vtk1, vtk2 : string
        The filenames of the vtk files to compare

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

    # construct element tree from vtk file
    root1 = ET.fromstring(open(vtk1).read())
    root2 = ET.fromstring(open(vtk2).read())

    # sort the vtk file in case nodes appear in different positions
    # e.g. because of minor changes in the output code
    sortedroot1 = sort_vtk(root1)
    sortedroot2 = sort_vtk(root2)

    if verbose:
        print("Comparing {} and {}".format(vtk1, vtk2))
        print("... with a maximum relative error of {} and a maximum absolute error of {}*max_abs_parameter_value.".format(relative, absolute))

    # sort the vtk file so that the comparison is independent of the
    # index numbering (coming e.g. from different grid managers)
    sortedroot1, sortedroot2 = sort_vtk_by_coordinates(sortedroot1, sortedroot2, verbose)

    # do the fuzzy compare
    if is_fuzzy_equal_node(sortedroot1, sortedroot2, absolute, relative, zeroValueThreshold, verbose):
        return 0
    else:
        return 1


# fuzzy compare of VTK nodes
def is_fuzzy_equal_node(node1, node2, absolute, relative, zeroValueThreshold, verbose):

    is_equal = True
    for node1child, node2child in zip(node1.iter(), node2.iter()):
        if node1.tag != node2.tag:
            if verbose:
                print('The name of the node differs in: {} and {}'.format(node1.tag, node2.tag))
                is_equal = False
            else:
                return False
        if list(node1.attrib.items()) != list(node2.attrib.items()):
            if verbose:
                print('Attributes differ in node: {}'.format(node1.tag))
                is_equal = False
            else:
                return False
        if len(list(node1.iter())) != len(list(node2.iter())):
            if verbose:
                print('Number of children differs in node: {}'.format(node1.tag))
                is_equal = False
            else:
                return False
        if node1child.text or node2child.text:
            if not is_fuzzy_equal_text(node1child.text, node2child.text,
                                       node1child.attrib["Name"],
                                       int(node1child.attrib["NumberOfComponents"]),
                                       absolute, relative, zeroValueThreshold, verbose):
                if node1child.attrib["Name"] == node2child.attrib["Name"]:
                    if verbose:
                        is_equal = False
                    else:
                        return False
                else:
                    if verbose:
                        print('Comparing different parameters: {} and {}'.format(node1child.attrib["Name"], node2child.attrib["Name"]))
                        is_equal = False
                    else:
                        return False
    return is_equal


# fuzzy compare of text (in the xml sense) consisting of whitespace separated numbers
def is_fuzzy_equal_text(text1, text2, parameter, numComp, absolute, relative, zeroValueThreshold, verbose):
    list1 = text1.split()
    list2 = text2.split()
    # difference only in whitespace?
    if (list1 == list2):
        return True
    # compare number by number
    is_equal = True

    # first split the list into compononents
    lists1 = []
    lists2 = []
    parameters = []
    for i in range(0, numComp):
        lists1.append(list1[i::numComp])
        lists2.append(list2[i::numComp])
        if numComp > 1:
            parameters.append("{}_{}".format(parameter, i))
            # if zero threshold was set for all components one component inherits it from the parameter
            if parameter in zeroValueThreshold:
                zeroValueThreshold["{}_{}".format(parameter, i)] = zeroValueThreshold[parameter]
        else:
            parameters.append(parameter)

    for list1, list2, parameter in zip(lists1, lists2, parameters):
        # for verbose output
        max_relative_difference = 0.0
        message = ''

        # see inspiration, explanations in
        # https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
        # get the order of magnitude of the parameter by calculating the max
        floatList1 = [float(i) for i in list1]
        floatList2 = [float(i) for i in list2]

        # check for nan and inf
        for number1, number2 in zip(floatList1, floatList2):
            if math.isnan(number1) or math.isnan(number2):
                print('Parameter {} contains NaN!'.format(parameter))
                return False
            if math.isinf(number1) or math.isinf(number2):
                print('Parameter {} contains inf!'.format(parameter))
                return False

        # manipulate the data set for the sake of sensible comparison
        # if the parameter is listed in the zeroThreshold dictionary replace all float under threshold with zero.
        # only replace them with zero if the parameters in both lists are under the threshold. Otherwise we
        # compare a non-zero value with 0 later.
        if parameter in zeroValueThreshold:
            floatList1 = [0.0 if abs(i) < float(zeroValueThreshold[parameter]) and abs(j) < float(zeroValueThreshold[parameter])
                          else i for i, j in zip(floatList1, floatList2)]
            floatList2 = [0.0 if abs(i) < float(zeroValueThreshold[parameter]) and abs(j) < float(zeroValueThreshold[parameter])
                          else j for i, j in zip(floatList1, floatList2)]

        absFloatList1 = [abs(i) for i in floatList1]
        absFloatList2 = [abs(i) for i in floatList2]

        magnitude = max(max(absFloatList1), max(absFloatList2))
        minimal = min(min(absFloatList1), min(absFloatList2))

        for number1, number2 in zip(floatList1, floatList2):
            diff = abs(number1 - number2)
            largernumber = max(abs(number1), abs(number2))

            # If the absolute criterion is satisfied we consider the numbers equal...
            # scale the absolute tolerance with the magnitude of the parameter
            if diff <= absolute * magnitude:
                continue

            # ...if not check the relative criterion
            if diff <= largernumber * relative:
                continue
            else:
                # the numbers are not equal
                if verbose:
                    is_equal = False
                    if largernumber != 0.0:
                        if diff / largernumber > max_relative_difference:
                            max_relative_difference = diff / largernumber
                            message = 'Difference is too large: {:.2%} -> between: {} and {}'.format(max_relative_difference, number1, number2)
                else:
                    return False

        if verbose and max_relative_difference != 0.0:
            print('\nData differs in parameter: {}'.format(parameter))
            print(message)
            print('Info for {}: max_abs_parameter_value={} and min_abs_parameter_value={}.'.format(parameter, magnitude, minimal))
            if parameter in zeroValueThreshold:
                print('For parameter {} a zero value threshold of {} was given.'.format(parameter, zeroValueThreshold[parameter]))

    return is_equal


def sort_by_name(elem):
    name = elem.get('Name')
    if name:
        try:
            return str(name)
        except ValueError:
            return ''
    return ''


# sorts attributes of an item and returns a sorted item
def sort_attributes(item, sorteditem):
    attrkeys = sorted(item.keys())
    for key in attrkeys:
        sorteditem.set(key, item.get(key))


def sort_elements(items, newroot):
    items = sorted(items, key=sort_by_name)
    items = sorted(items, key=attrgetter('tag'))

    # Once sorted, we sort each of the items
    for item in items:
        # Create a new item to represent the sorted version
        # of the next item, and copy the tag name and contents
        newitem = ET.Element(item.tag)
        if item.text and item.text.isspace() == False:
            newitem.text = item.text

        # Copy the attributes (sorted by key) to the new item
        sort_attributes(item, newitem)

        # Copy the children of item (sorted) to the new item
        sort_elements(list(item), newitem)

        # Append this sorted item to the sorted root
        newroot.append(newitem)


# has to sort all Cell and Point Data after the attribute "Name"!
def sort_vtk(root):
    if(root.tag != "VTKFile"):
        print('Format is not a VTKFile. Sorting will most likely fail!')
    # create a new root for the sorted tree
    newroot = ET.Element(root.tag)
    # create the sorted copy
    # (after the idea of Dale Lane's xmldiff.py)
    sort_attributes(root, newroot)
    sort_elements(list(root), newroot)
    # return the sorted element tree
    return newroot


# sorts the data by point coordinates so that it is independent of index numbering
def sort_vtk_by_coordinates(root1, root2, verbose):
    if not is_fuzzy_equal_node(root1.find(".//Points/DataArray"), root2.find(".//Points/DataArray"), absolute=1e-2, relative=1.5e-7, zeroValueThreshold=dict(), verbose=False):
        if verbose:
            print("Sorting vtu by coordinates...")
        for root in [root1, root2]:
            # parse all DataArrays in a dictionary
            pointDataArrays = []
            cellDataArrays = []
            dataArrays = {}
            numberOfComponents = {}
            for dataArray in root.findall(".//PointData/DataArray"):
                pointDataArrays.append(dataArray.attrib["Name"])
            for dataArray in root.findall(".//CellData/DataArray"):
                cellDataArrays.append(dataArray.attrib["Name"])
            for dataArray in root.findall(".//DataArray"):
                dataArrays[dataArray.attrib["Name"]] = dataArray.text
                numberOfComponents[dataArray.attrib["Name"]] = dataArray.attrib["NumberOfComponents"]

            vertexArray = []
            coords = dataArrays["Coordinates"].split()
            # group the coordinates into coordinate tuples
            dim = int(numberOfComponents["Coordinates"])
            for i in range(len(coords) // dim):
                vertexArray.append([float(c) for c in coords[i * dim: i * dim + dim]])

            # group the cells into vertex index tuples
            cellArray = []
            offsets = dataArrays["offsets"].split()
            connectivity = dataArrays["connectivity"].split()
            vertex = 0
            for cellIdx, offset in enumerate(offsets):
                cellArray.append([])
                for v in range(vertex, int(offset)):
                    cellArray[cellIdx].append(int(connectivity[v]))
                    vertex += 1

            # for non-conforming output vertices can have the same coordinates and also
            # different indices / sorting so we need another criterium to sort.
            # we use the largest cell midpoint coordinate vector the vertex is connected to
            largestCellMidPointForVertex = [[0, 0, 0]]*len(vertexArray)
            for cellIdx, cell in enumerate(cellArray):
                # compute cell midpoint
                coords = [vertexArray[i] for i in cell]
                midpoint = [i/float(len(coords)) for i in [sum(coord) for coord in zip(*coords)]]
                for vertexIndex in cell:
                    largestCellMidPointForVertex[vertexIndex] = max(largestCellMidPointForVertex[vertexIndex], midpoint)

            # obtain a vertex index map
            vMap = []
            for idx, coords in enumerate(vertexArray):
                vMap.append((idx, coords))

            vertexIndexMap = [0]*len(vMap)
            vertexIndexMapInverse = [0]*len(vMap)
            # first sort by coordinates, if the same by largestCellMidPointForVertex
            for idxNew, idxOld in enumerate(sorted(vMap, key=lambda x: (x[1], largestCellMidPointForVertex[x[0]]))):
                vertexIndexMap[idxOld[0]] = idxNew
                vertexIndexMapInverse[idxNew] = idxOld[0]

            # replace all vertex indices in the cellArray by the new indices
            for i, cell in enumerate(cellArray):
                for j, vertexIndex in enumerate(cell):
                    cellArray[i][j] = vertexIndexMap[vertexIndex]

            # sort the indices of all cells in case the two grids have different orientation of the elements
            # this doesn't necessarily have to be a valid vtk order, it's just temporary for a sorted literal comparison
            # of the files
            for i, cell in enumerate(cellArray):
                cellArray[i] = sorted(cell)

            # sort all data arrays
            for name, text in list(dataArrays.items()):
                # split the text
                items = text.split()
                # convert if vector
                num = int(numberOfComponents[name])
                newitems = []
                for i in range(len(items) // num):
                    newitems.append([i for i in items[i * num: i * num + num]])
                items = newitems
                # sort the items: we have either vertex or cell data
                if name in pointDataArrays:
                    sortedItems = [items[i] for i in vertexIndexMapInverse]
                elif name in cellDataArrays or name == "types":
                    sortedItems = [j for (i, j) in sorted(zip(cellArray, items), key=itemgetter(0))]
                elif name == "offsets":
                    sortedItems = []
                    counter = 0
                    for cell in sorted(cellArray):
                        counter += len(cell)
                        sortedItems.append([str(counter)])
                elif name == "Coordinates":
                    sortedItems = [vertexArray[i] for i in vertexIndexMapInverse]
                elif name == "connectivity":
                    sortedItems = sorted(cellArray)

                # convert the sorted arrays to a xml text
                dataArrays[name] = ""
                for i in sortedItems:
                    for j in i:
                        dataArrays[name] += str(j) + " "

            # do the replacement in the actual elements
            for dataArray in root.findall(".//DataArray"):
                dataArray.text = dataArrays[dataArray.attrib["Name"]]

    return (root1, root2)


# main program if called as script return appropriate error codes
if __name__ == "__main__":
    # handle arguments and print help message
    parser = argparse.ArgumentParser(description='Fuzzy compare of two VTK\
        (Visualization Toolkit) files. The files are accepted if for every\
        value the difference is below the absolute error or below the\
        relative error or below both.')
    parser.add_argument('vtk_file_1', type=str, help='first file to compare')
    parser.add_argument('vtk_file_2', type=str, help='second file to compare')
    parser.add_argument('-r', '--relative', type=float, default=1e-2, help='maximum relative error (default=1e-2)')
    parser.add_argument('-a', '--absolute', type=float, default=1.5e-7, help='maximum absolute error (default=1.5e-7)')
    parser.add_argument('-z', '--zeroThreshold', type=json.loads, default='{}', help='Thresholds for treating numbers as zero for a parameter as a python dict e.g. {"vel":1e-7,"delP":1.0}')
    parser.add_argument('-v', '--verbose', type=bool, default=True, help='verbosity of the script')
    args = vars(parser.parse_args())

    sys.exit(compare_vtk(args["vtk_file_1"], args["vtk_file_2"], args["absolute"], args["relative"], args["zeroThreshold"], args["verbose"]))
