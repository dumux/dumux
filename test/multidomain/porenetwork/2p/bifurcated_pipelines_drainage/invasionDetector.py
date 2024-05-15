import glob, os
import re
import subprocess
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from generate_network_analyze import generate_network_with_random_throat_radii
from generate_network_analyze import pcEntry_invasion_path

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
    invaded_throat_indices = sorted(invaded_throat_indices, key=lambda x: pcEntry_invasion_path(x))
    return invaded_throat_indices

def get_pc_for_invading_pores(invading_path, number_invaded, vtp_files):
    invasion_sequence = np.array_split(invading_path, np.cumsum(number_invaded))
    invasion_pc_ij = []
    reference_pc_ij = []
    reader = vtk.vtkXMLPolyDataReader()
    for invaded_throats, vtp in zip(invasion_sequence, vtp_files[1:]):
        reader.SetFileName(vtp)
        reader.Update()
        polydata = reader.GetOutput()
        pc_ij_array = polydata.GetCellData().GetArray('pcij')
        pc_ij = vtk_to_numpy(pc_ij_array)
        for idx in invaded_throats:
            if (idx == invaded_throats[-1]):
                invasion_pc_ij.append(pc_ij[idx])
                reference_pc_ij.append(pcEntry_invasion_path(idx))
    return [reference_pc_ij, invasion_pc_ij]

if __name__ == "__main__":
    generate_network_and_analyze = 'python3 generate_network_analyze.py' #cmd to generate network
    simulation_runs_for_each_time_step_size = 500

    testName = ['validate_preferential_path', 'validate_preferential_path_reg', 'validate_preferential_path_md']
    maxTimeStep =  [16, 8, 4, 2, 1] # time steps which are going to be used
    regularizationDelta = [0.4, 0.1, 1e-3] # we test few different regularization Delta

    prediction_rate = {}
    avg_l2error_pcentry = {}


    for test in testName:                   # initialize dict with test name as key
        if (test is not testName[1]):
            prediction_rate[test] = []
            avg_l2error_pcentry[test] = []
        else:
            for regDelta in regularizationDelta:
                prediction_rate[test + str(regDelta)] = []
                avg_l2error_pcentry[test + str(regDelta)] = []

    for maxDtIdx, maxDt in enumerate(maxTimeStep):
        counter_correct_prediction = {}   # count how many times prediction is correct
        l2errorPcCollector = {}


        for test in testName:
            if (test is not testName[1]):  # FI-N or FI-Theta
                l2errorPcCollector[test] = []
                counter_correct_prediction[test] = 0
            else:
                for regDelta in regularizationDelta:
                        l2errorPcCollector[test + str(regDelta)] =[]
                        counter_correct_prediction[test + str(regDelta)] = 0


        for i in range(simulation_runs_for_each_time_step_size):
            generate_network_with_random_throat_radii(i)                    # generate a network
            analyticalSolution = np.load("analytical_solution.npz")         # load anaylitcal solution
            analytical_path = analyticalSolution['analytical_path']         # analytical invasion order
            analytical_pcEntry = analyticalSolution['analytical_pcEntry']   # analytical pc entry for invasion

            for test in testName:
                if (test is not testName[1]):  # FI-N or FI-Theta
                    if (test is testName[0]):
                        param = 'params.input'
                    elif (test is testName[2]):
                        param = 'params_md.input'
                    delete_all_vtp_files()
                    subprocess.run(['./' + test] + [param]
                                   + ['-TimeLoop.DtInitial', str(maxDt)]
                                   + ['-TimeLoop.TEnd', str(1500)]
                                   + ['-TimeLoop.MaxTimeStepSize', str(maxDt)]
                                   + ['-Grid.ThroatCrossSectionShape', 'Circle']
                                   + ['-Problem.Name', str(test) + '_run_' + str(i)]
                                   + ['-Pnm.Problem.Name', str(test) + '_run_' + str(i)]
                                   + ['-Problem.VtpOutputFrequency', str(0)])
                    numerical_invasion_order = np.array([])
                    vtp_file_lists = find_all_vtp_files_withPrefix(test + '_run_' + str(i))
                    vtp_file_lists.sort(key=lambda x: extract_number(x))
                    invasionStates = []
                    number_invaded_throats = []
                    for vtp in vtp_file_lists:
                        invasionStates.append(get_invasion_status_from_vtp(vtp))
                    for previous, current in zip(invasionStates, invasionStates[1:]):
                        invaded_throat_indices = get_new_invaded_throat_indices(previous, current)
                        number_invaded_throats.append(len(invaded_throat_indices))
                        numerical_invasion_order = np.append(numerical_invasion_order, invaded_throat_indices)
                    numerical_invasion_order = numerical_invasion_order.astype(int)
                    correctpredition = np.array_equal(numerical_invasion_order, analytical_path)
                    print(analytical_path)
                    print(analytical_pcEntry)
                    print(numerical_invasion_order)
                    if (correctpredition):
                        counter_correct_prediction[test] = counter_correct_prediction[test] + 1
                    [reference_pcEntry, numerical_pcEntry] = get_pc_for_invading_pores(numerical_invasion_order, number_invaded_throats, vtp_file_lists)
                    print(numerical_pcEntry)
                    print(correctpredition)
                    if (correctpredition):
                        l2errorPc = np.sum(np.power((np.array(numerical_pcEntry)-np.array(reference_pcEntry)),2))/len(reference_pcEntry)
                        l2errorPc = np.sqrt(l2errorPc)
                        l2errorPcCollector[test].append(l2errorPc)
                    # Write all important info of a simulation run into a txt file
                    with open('results_' + str(test) + '_run_' + str(i) + '_maxDt_' +  str(maxDt) + '.txt', 'w') as f:
                        f.write("** The analytical throat invasion order is: **\n")
                        f.write(str(analytical_path) + '\n')
                        f.write("** The analytical entry pressure for invaded throats are: **\n")
                        f.write(str(analytical_pcEntry) + '\n')
                        f.write("** The throat invasion order numerically calculated is: **\n")
                        f.write(str(numerical_invasion_order) + '\n')
                        f.write("** The pc_ij at throat for the time step where invasion really captured: **\n")
                        f.write(str(numerical_pcEntry) + '\n')
                        f.write("** The reference pc_ij at throat compared to the numerical pc_ij: **\n")
                        f.write(str(reference_pcEntry) + '\n')
                        if (np.array_equal(numerical_invasion_order, analytical_path)):
                            f.write("The invasion order prediction for this simulation run is correct.\n")
                            f.write("The L2 error of pc_ij and pc_entry is : " + str(l2errorPc))
                        else:
                            f.write("The invasion order prediction for this simulation run is wrong.\n")

                elif (test is testName[1]): # FI-R, loop for different Delta
                    for regDelta in regularizationDelta:
                        delete_all_vtp_files()
                        subprocess.run(['./' + test]
                                       + ['-TimeLoop.DtInitial', str(maxDt)]
                                       + ['-TimeLoop.TEnd', str(1500)]
                                       + ['-TimeLoop.MaxTimeStepSize', str(maxDt)]
                                       + ['-Grid.ThroatCrossSectionShape', 'Circle']
                                       + ['-Pnm.Problem.Name', str(test) + '_regdelta_' + str(regDelta) + '_run_' + str(i)]
                                       + ['-Problem.RegularizationDelta', str(regDelta/(2**maxDtIdx))]
                                       + ['-Problem.VtpOutputFrequency', str(0)])
                        numerical_invasion_order = np.array([])
                        vtp_file_lists = find_all_vtp_files_withPrefix(str(test) + '_regdelta_' + str(regDelta) + '_run_' + str(i))
                        vtp_file_lists.sort(key=lambda x: extract_number(x))
                        invasionStates = []
                        number_invaded_throats = []
                        for vtp in vtp_file_lists:
                            invasionStates.append(get_invasion_status_from_vtp(vtp))
                        for previous, current in zip(invasionStates, invasionStates[1:]):
                            invaded_throat_indices = get_new_invaded_throat_indices(previous, current)
                            number_invaded_throats.append(len(invaded_throat_indices))
                            numerical_invasion_order = np.append(numerical_invasion_order, invaded_throat_indices)
                        numerical_invasion_order = numerical_invasion_order.astype(int)
                        correctpredition = np.array_equal(numerical_invasion_order, analytical_path)
                        if (correctpredition):
                            counter_correct_prediction[test + str(regDelta)] = counter_correct_prediction[test + str(regDelta)] + 1
                        [reference_pcEntry, numerical_pcEntry] = get_pc_for_invading_pores(numerical_invasion_order, number_invaded_throats, vtp_file_lists)
                        if (correctpredition):
                            l2errorPc = np.sum(np.power((np.array(numerical_pcEntry)-np.array(reference_pcEntry)),2))/len(reference_pcEntry)
                            l2errorPc = np.sqrt(l2errorPc)
                            l2errorPcCollector[test + str(regDelta)].append(l2errorPc)
                        # Write all important info of a simulation run into a txt file
                        with open('results_' + str(test) + '_regdelta_' + str(regDelta) + '_run_' + str(i) + '_maxDt_' +  str(maxDt) + '.txt', 'w') as f:
                            f.write("** The analytical throat invasion order is: **\n")
                            f.write(str(analytical_path) + '\n')
                            f.write("** The analytical entry pressure for invaded throats are: **\n")
                            f.write(str(analytical_pcEntry) + '\n')
                            f.write("** The throat invasion order numerically calculated is: **\n")
                            f.write(str(numerical_invasion_order) + '\n')
                            f.write("** The pc_ij at throat for the time step where invasion really captured: **\n")
                            f.write(str(numerical_pcEntry) + '\n')
                            f.write("** The reference pc_ij at throat compared to the numerical pc_ij: **\n")
                            f.write(str(reference_pcEntry) + '\n')

                            if (np.array_equal(numerical_invasion_order, analytical_path)):
                                f.write("The invasion order prediction for this simulation run is correct.\n")
                                f.write("The L2 error of pc_ij and pc_entry is : " + str(l2errorPc))
                            else:
                                f.write("The invasion order prediction for this simulation run is wrong.\n")

        for test in testName:
            if (test is not testName[1]):
                prediction_rate[test].append(counter_correct_prediction[test] / simulation_runs_for_each_time_step_size)
                avg_l2error_pcentry[test].append(np.average(l2errorPcCollector[test]))
            else:
                for delta in regularizationDelta:
                    print("counter correct prediction for", delta, "is", counter_correct_prediction[test + str(delta)])
                    prediction_rate[test + str(delta)].append(counter_correct_prediction[test + str(delta)] / simulation_runs_for_each_time_step_size)
                    avg_l2error_pcentry[test + str(delta)].append(np.average(l2errorPcCollector[test + str(delta)]))



fig, ax = plt.subplots(dpi=300, ncols=2, nrows=1, figsize=(6, 3))

methods = ['FI-N', 'FI-R', 'FI-$\Theta$']
for method, test in zip(methods, testName):
    if (test is not testName[1]):
        ax[0].plot(maxTimeStep, prediction_rate[test], label=method)
        ax[1].plot(maxTimeStep, avg_l2error_pcentry[test], label = method)
    elif (test is testName[1]):
        for delta in regularizationDelta:
            ax[0].plot(maxTimeStep, prediction_rate[test + str(delta)], label = method + ' $\delta$ = ' + str(delta))
            ax[1].plot(maxTimeStep, avg_l2error_pcentry[test + str(delta)], label = method + ' $\delta$ = ' + str(delta))
ax[0].set_xlabel("Maximum time step size")
ax[0].set_ylabel("Accuracy of prediction [-]")
ax[1].set_xlabel("Maximum time step size")
ax[1].set_ylabel("$E_{p_{c,e}}$ [-]")
ax[1].set_yscale('log')
ax[0].legend()
ax[1].legend()
plt.subplots_adjust(wspace=0, hspace=0)
fig.tight_layout(rect=[0.03, 0.07, 1, 0.9], pad=0.4, w_pad=2.0, h_pad=1.0)
plt.show()
plt.savefig("Random_test.pdf", dpi=900)
