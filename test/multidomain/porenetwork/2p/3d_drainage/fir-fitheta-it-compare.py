import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

def simulation_succeed(tEnd):
    last_time_step = np.genfromtxt("logfile_pnm_3d_md.txt", usecols = 0).T[-1]
    return abs(last_time_step - tEnd) < 1e-3

if __name__ == "__main__":
    simulation_runs = 20
    tEnd = 30.0
    no_scaling_it = [] # list to store iteration numbers not using theta-scaling
    scaling_it = [] # list to store iteration numbers using theta-scaling
    fi_r_it = []
    for i in range(simulation_runs):
        subprocess.call(['python3', 'generate_lattice_grid.py'])
        print("Now we test grid ", i, "...")
        print("First we test FI-theta without any scaling")
        subprocess.run(  ['./pnm_3d_drainage_md', 'params.input']
                       + ['-Grid.File', 'lattice_network.dgf']
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-Problem.NonWettingMassFlux', str(5e-8)]
                       + ['-Newton.MaxSteps', str(10)]
                       + ['-Newton.TargetSteps', str(4)]
                       + ['-Newton.MaxTimeStepDivisions', str(20)]
                       + ['-Pnm.Problem.Name', 'pnm_3d_md']
                       + ['-Constraint.Problem.Name', 'pnm_3d_md_constraint']
                       + ['-Constraint.Problem.ThetaScalingFactor', str(1)])
        if simulation_succeed(tEnd):
            print('simulation succeeded for grid ', i, ' without using scaling')
            no_scaling_it.append(np.genfromtxt('NewtonLog.txt'))
        else:
            print('simulation failed for grid ', i, ' without using scaling')
            no_scaling_it.append(0)

        print("Then we test FI-theta with theta-scaling: ")
        subprocess.run(  ['./pnm_3d_drainage_md', 'params.input']
                       + ['-Grid.File', 'lattice_network.dgf']
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-Problem.NonWettingMassFlux', str(5e-8)]
                       + ['-Newton.MaxSteps', str(10)]
                       + ['-Newton.TargeTEndtSteps', str(4)]
                       + ['-Newton.MaxTimeStepDivisions', str(20)]
                       + ['-Pnm.Problem.Name', 'pnm_3d_md']
                       + ['-Constraint.Problem.Name', 'pnm_3d_md_constraint']
                       + ['-Constraint.Problem.ThetaScalingFactor', str(1e-1)])
        if simulation_succeed(tEnd):
            print('simulation succeeded for grid ', i, ' using scaling')
            scaling_it.append(np.genfromtxt('NewtonLog.txt'))
        else:
            print('simulation failed for grid ', i, ' using scaling')
            scaling_it.append(0)

        print("Now we test fi-r: ")
        subprocess.run(  ['./pnm_3d_drainage_reg', 'params.input']
                       + ['-Grid.File', 'lattice_network.dgf']
                       + ['-Problem.NonWettingMassFlux', str(5e-8)]
                       + ['-TimeLoop.TEnd', str(tEnd)]
                       + ['-Newton.MaxSteps', str(10)]
                       + ['-Newton.TargetSteps', str(4)]
                       + ['-Newton.MaxTimeStepDivisions', str(20)]
                       + ["-Problem.RegularizationDelta", str(0.01)])
        if simulation_succeed(tEnd):
            print('simulation succeeded for grid ', i, ' using fi-r')
            fi_r_it.append(np.genfromtxt('NewtonLog.txt'))
        else:
            print('simulation failed for grid ', i, ' using fi-r')
            fi_r_it.append(0)

    plt.plot(np.arange(simulation_runs), no_scaling_it, marker = 'o', markersize = 10, label ="without $\Theta$ scaling")
    plt.plot(np.arange(simulation_runs), scaling_it, marker = 'x', markersize = 10, label = "$\Theta$-scaling, C=0.1")
    plt.plot(np.arange(simulation_runs), fi_r_it, marker = '>', markersize = 10, label = "FI-R, $\tilde{\theta}$ = 0.01")
    plt.xlabel('simulation runs')
    plt.ylabel('iteration numbers')
    plt.legend()
    plt.show()
