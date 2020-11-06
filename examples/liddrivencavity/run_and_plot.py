#!/usr/bin/env python3
import sys
if sys.version_info[0] < 3:
    sys.exit("Python 3 required to run this script. Exit.")

#####################################################
##### simulation ####################################
#####################################################
import subprocess
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='plot script for the lid-driven cavity example')
parser.add_argument('-s', '--skipsim', required=False, action='store_true',
                    help='use this flag to skip the simulation run and directly plot already existing data')
parser.add_argument('-n', '--noplotwindow', required=False, action='store_true',
                    help='use this flag to suppress the plot window popping up')
args = vars(parser.parse_args())

reynolds = [1, 1000]
x = {}
y = {}
vx = {}
vy = {}
for re in reynolds:
    if not args['skipsim']:
        subprocess.run(['make', 'example_ff_liddrivencavity'], check=True)
        subprocess.run(['./example_ff_liddrivencavity', 'params_re' + str(re) + '.input'], check=True)

    y[str(re)], vx[str(re)] = np.genfromtxt('example_ff_liddrivencavity_re' + str(re) + '_vx' + '.log', skip_header= True).T
    x[str(re)], vy[str(re)] = np.genfromtxt('example_ff_liddrivencavity_re' + str(re) + '_vy' + '.log', skip_header= True).T

####################################################
#### reference #####################################
####################################################
ghiavx, ghiay = np.genfromtxt("./reference_data/ghia_x.csv", delimiter=';').T
ghiax, ghiavy = np.genfromtxt("./reference_data/ghia_y.csv", delimiter=';').T

jurjevicnumvx, jurjevicnumy = np.genfromtxt("./reference_data/v_x_num.csv", delimiter=';').T
jurjevicnumy = jurjevicnumy/2.0 + 0.5
jurjevicnumx, jurjevicnumvy = np.genfromtxt("./reference_data/v_y_num.csv", delimiter=';').T
jurjevicnumx = jurjevicnumx/2.0 + 0.5
jurjevicexpvx, jurjevicexpy = np.genfromtxt("./reference_data/v_x_exp.csv", delimiter=';').T
jurjevicexpy = jurjevicexpy/2.0 + 0.5
jurjevicexpx, jurjevicexpvy = np.genfromtxt("./reference_data/v_y_exp.csv", delimiter=';').T
jurjevicexpx = jurjevicexpx/2.0 + 0.5

####################################################
#### plotting ######################################
####################################################
try:
    import matplotlib
    import matplotlib.pyplot as plt

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (9,4))
    ax1.plot(vx['1'], y['1'] , color='black', label=u"DuMu$^\mathrm{x}$",linewidth=2)
    ax1.plot(jurjevicnumvx, jurjevicnumy, '--', markerfacecolor='white', color='black', label=u"R.Jurjevic et al., num")
    ax1.plot(jurjevicexpvx, jurjevicexpy, 'o', markerfacecolor='white', color='black', label=u"R.Jurjevic et al., exp")
    ax1.set_xlabel(r"$v_x$[m/s]")
    ax1.set_ylabel(u"y [m]")

    ax2.plot(x['1'], vy['1'],  color='black', label=u"DuMu$^\mathrm{x}$",linewidth=2)
    ax2.plot(jurjevicnumx, jurjevicnumvy, '--', markerfacecolor='white', color='black', label=u"R.Jurjevic, num")
    ax2.plot(jurjevicexpx, jurjevicexpvy, 'o', markerfacecolor='white', color='black', label=u"R.Jurjevic, exp")
    ax2.set_xlabel(u"x [m]")
    ax2.set_ylabel(r"$v_y$[m/s]",labelpad=1)
    ax2.set_xlabel(u"x [m]")
    ax2.set_ylabel(r"$v_y$[m/s]",labelpad=1)

    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(0.51, 1.0), ncol=3, labelspacing=0.)

    ax3.plot(vx['1000'], y['1000'], color='black', label=u"DuMu$^\mathrm{x}$",linewidth=2)
    ax3.plot(ghiavx, ghiay, 'o',  markerfacecolor='white', color='black', label=u"Ghia et al.")
    ax3.set_xlabel(r"$v_x$[m/s]")
    ax3.set_ylabel(u"y [m]")

    ax4.plot(x['1000'], vy['1000'], color='black', label=u"DuMu$^\mathrm{x}$",linewidth=2)
    ax4.plot(ghiax, ghiavy, 'o',  markerfacecolor='white', color='black', label=u"Ghia et al.")
    ax4.set_xlabel(u"x [m]")
    ax4.set_ylabel(r"$v_y$[m/s]",labelpad=1)

    handles, labels = ax4.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(0.92, 1.0), ncol =2, labelspacing=0.)
    fig.tight_layout(rect=[0.03, 0.07, 1, 0.9], pad=0.4, w_pad=2.0, h_pad=1.0)

    plt.savefig("lidverification.png", dpi= 300)
    if not args['noplotwindow']: plt.show()

except ImportError:
    print("Skipping plot: matplotlib has not been found.")
