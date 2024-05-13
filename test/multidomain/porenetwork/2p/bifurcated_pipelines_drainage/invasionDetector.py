import glob, os
import re
import subprocess
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def find_all_vtp_files_withPrefix(prefix):
    vtp_file_lists = []
    for file in glob.glob(prefix + "*.vtp"):
        vtp_file_lists.append(file)
    return vtp_file_lists

def extract_number(f):
    s = re.findall("(\d+).vtp",f)
    return (int(s[0]) if s else -1,f)

def find_last_vtp_file(vtp_file_lists):
    return max(vtp_file_lists, key=extract_number)

def delete_all_vtp_files():
    for item in os.listdir():
        if item.endswith(".vtp"):
            os.remove(item)

def get_time_steps_from_pvd(pvd_file):
    time_steps = []
    tree = ET.parse(pvd_file)
    root = tree.getroot()
    for collection in root.findall('.//DataSet'):
        timestep_attr = collection.attrib.get('timestep')
        if timestep_attr is not None:
            time_step = float(timestep_attr)
            time_steps.append(time_step)
    return time_steps

def get_invasion_status_from_vtp(vtp_file_name):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_file_name)
    reader.Update()
    polydata = reader.GetOutput()
    invaded_vtk = polydata.GetCellData().GetArray('Invaded')
    invaded = vtk_to_numpy(invaded_vtk)
    return invaded

def get_new_invaded_throat_indices(previous, current):
    diff = (~np.equal(previous, current)).astype(int)
    invaded_throat_indices = np.flatnonzero(diff)
    return invaded_throat_indices

def get_pc_for_invading_pores(invading_path, vtp_files):
    invasion_pc_ij = []
    reader = vtk.vtkXMLPolyDataReader()
    for idx, vtp in zip(invading_path, vtp_files[1:]):
        reader.SetFileName(vtp)
        reader.Update()
        polydata = reader.GetOutput()
        pc_ij_array = polydata.GetCellData().GetArray('pcij')
        pc_ij = vtk_to_numpy(pc_ij_array)
        invasion_pc_ij.append(pc_ij[idx])
    return invasion_pc_ij

if __name__ == "__main__":
    analyticalSolution = np.load("analytical_solution.npz")
    analytical_path = analyticalSolution['analytical_path']
    print("The analyitcal invasion path is: ", analytical_path)
    analytical_pcEntry = analyticalSolution['analytical_pcEntry']
    print("The cooresponding Pc entry is : ", analytical_pcEntry)
    numerical_invasion_order = np.array([])
    vtp_file_lists = find_all_vtp_files_withPrefix("validate_preferential_path")
    vtp_file_lists.sort(key=lambda x: extract_number(x))
    invasionStates = []
    for vtp in vtp_file_lists:
        invasionStates.append(get_invasion_status_from_vtp(vtp))
    for previous, current in zip(invasionStates, invasionStates[1:]):
        invaded_throat_indices = get_new_invaded_throat_indices(previous, current)
        numerical_invasion_order = np.append(numerical_invasion_order, invaded_throat_indices)

    numerical_invasion_order = numerical_invasion_order.astype(int)
    print(numerical_invasion_order)
    if (numerical_invasion_order == numerical_invasion_order).all():
        print("The prediction is correct!")

    numerical_pcEntry = get_pc_for_invading_pores(numerical_invasion_order, vtp_file_lists)
    print("Numerical pc entry: ", numerical_pcEntry)
