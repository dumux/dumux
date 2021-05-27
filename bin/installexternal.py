#!/usr/bin/env python3

"""
install external stuff for dumux
"""
import os
import urllib.request
import tarfile
import sys
import subprocess
import shutil
import re
import argparse

class ChoicesAction(argparse._StoreAction):
    def __init__(self, **kwargs):
        super(ChoicesAction, self).__init__(**kwargs)
        if self.choices is None:
            self.choices = []
        self._choices_actions = []
    def add_choice(self, choice, help=''):
        self.choices.append(choice)
        choice_action = argparse.Action(option_strings=[], dest=choice, help=help)
        self._choices_actions.append(choice_action)
    def _get_subactions(self):
        return self._choices_actions

def show_message(message):
    print("*" * 120)
    print(message)
    print("")
    print("*" * 120)


if len(sys.argv) == 1:
    show_message('No options given. For more information run the following command: \n ./installexternal.py --help')
    sys.exit()

parser = argparse.ArgumentParser(prog='installexternal',
                                 usage='./installexternal.py [OPTIONS] PACKAGES',
                                 description='This script downloads extenstions for dumux and dune \
                                     and installs some External Libraries and Modules.')
parser.register('action', 'store_choice', ChoicesAction)
# Positional arguments
group = parser.add_argument_group(title='your choice of packages')
packages = group.add_argument('packages', nargs='+', metavar='PACKAGES',
                 action='store_choice')
packages.add_choice('dumux-extensions', help="Download dumux-course and dumux-lecture.")
packages.add_choice('dune-extensions', help="Download dune-uggrid, dune-alugrid, dune-foamgrid, \
                    dune-subgrid, dune-spgrid, dune-mmesh and dune-functions.")
packages.add_choice('optimization', help="Download and install glpk and nlopt.")
packages.add_choice('others', help="Download and install opm , metis and gstat.")
packages.add_choice('lecture', help="Download dumux-lecture.")
packages.add_choice('course', help="Download dumux-course.")
packages.add_choice('ug', help="Download dune-uggrid.")
packages.add_choice('alugrid', help="Download dune-alugrid.")
packages.add_choice('foamgrid', help="Download dune-foamgrid.")
packages.add_choice('subgrid', help="Download dune-subgrid.")
packages.add_choice('spgrid', help="Download dune-spgrid.")
packages.add_choice('mmesh', help="Download dune-mmesh.")
packages.add_choice('functions', help="Download dune-functions.")
packages.add_choice('glpk', help="Download and install glpk.")
packages.add_choice('nlopt', help="Download and install nlopt.")
packages.add_choice('opm', help="Download opm modules required for cornerpoint grids.")
packages.add_choice('metis', help="Download and install the METIS graph partitioner.")
packages.add_choice('gstat', help="Download and install gstat.")


# Optional arguments
options = parser.add_mutually_exclusive_group(required=False)
options.add_argument('--clean', action="store_true", default=False,
                     help='Delete all files for the given packages.')
options.add_argument('--download', action="store_true", default=False,
                     help='Only download the packages.')

parser.add_argument('--dune_branch', default="releases/2.7",
                     help='Dune branch to be checked out.')
parser.add_argument('--dumux_branch', default="releases/3.4",
                     help='Dumux branch to be checked out.')
parser.add_argument('--opm_branch', default="release/2020.10",
                     help='Opm branch to be checked out.')
parser.add_argument('--mmesh_branch', default="release/1.2",
                     help='Mmesh branch to be checked out.')

args = vars(parser.parse_args())

def run_command(command, currentdir='.'):
    with open(currentdir+"/installexternal.log", "a") as log:
        popen = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for line in popen.stdout:
            log.write(line)
            print(line, end='')
        for line in popen.stderr:
            log.write(line)
            print(line, end='')
        popen.stdout.close()
        popen.stderr.close()
        return_code = popen.wait()
        if return_code:
            print("\n")
            message = "\n    (Error) The command {} returned with non-zero exit code\n".format(command)
            message += "\n    If you can't fix the problem yourself consider reporting your issue\n"
            message += "    on the mailing list (dumux@listserv.uni-stuttgart.de) and attach the file 'installexternal.log'\n"
            show_message(message)
            sys.exit(1)

def git_clone(url, branch=None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    result = run_command(command=[*clone, url])

def install_external(args):
    dune_branch = args['dune_branch']
    dumux_branch = args['dumux_branch']
    opm_branch = args['opm_branch']
    mmesh_branch = args['mmesh_branch']
    packages = args['packages']
    cleanup = args['clean']
    download = args['download']

    final_message = []
    top_dir = os.getcwd()
    ext_dir =  top_dir + "/external"


    # Prepare a list of packages
    packages = []
    for pkg in args['packages']:
        if pkg in packagenames:
            packages.extend(packagenames[pkg])
        else:
            packages.extend([key for key in external_urls.keys() if pkg in key])
    args['packages'] = packages

    # print the list of packages to be downloaded/installed/removed
    print("The following package(s) will be {0}:\n".format('removed' if cleanup else 'downloaded'), ', '.join(args['packages']), "\n")

    # check Location For DuneModules
    if not os.path.isdir("dune-common"):
        show_message(
            "You have to call " + sys.argv[0] + " for " + sys.argv[1] + " from\n"
            "the same directory in which dune-common is located.\n"
            "You cannot install it in this folder.")
        return

    # clear the log file
    logdir = ext_dir if os.path.exists(ext_dir) else top_dir
    open(logdir+'/installexternal.log', 'w').close()

    for package in packages:
        os.chdir(top_dir)
        # Package name for final message
        final_message.append('[---'+package+'---]')

        # Set the directory: create ext_dir for external packages
        if not any([re.compile(p).match(package) for p in ['dumux','dune', 'opm']]):
            os.makedirs(ext_dir, exist_ok=True)
            os.chdir(ext_dir)

        # Set the branch
        if 'dumux' in package:
            branch = dumux_branch
        elif 'mmesh' in package:
            branch = mmesh_branch
        elif 'dune' in package:
            branch = dune_branch
        elif 'opm' in package:
            branch = opm_branch

        # Run the requested command
        if cleanup:
            if os.path.isfile(package + '.tar.gz'):
                os.remove(package + '.tar.gz')
            if os.path.exists(package):
                # Remove
                shutil.rmtree(package)

                # Save message to be shown at the end
                final_message.append("{} has been removed.".format(package))
            else:
                # Save message to be shown at the end
                final_message.append("The folder {} does not exist.".format(package))
            continue

        else:
            # Check if tarball
            tarball = external_urls[package].endswith('tar.gz')

            if not os.path.exists(package):

                if tarball:

                    # Download the tarfile
                    filedata = urllib.request.urlopen(external_urls[package])
                    datatowrite = filedata.read()

                    with open(ext_dir + "/" + package +".tar.gz", 'wb') as f:
                        f.write(datatowrite)
                    # Save message to be shown at the end
                    final_message.append("{} has been sucessfully downloaded.".format(package))

                    # Start Installation if the flag download is set to false.
                    if not download:
                        # Extract
                        tf = tarfile.open(package+".tar.gz")
                        tf.extractall()

                        # Rename
                        shutil.move(os.path.commonprefix(tf.getnames()), package)

                        # Start the configuration
                        os.chdir(ext_dir + "/"+package)
                        if package == 'gstat':
                            with open('configure', 'r+') as f:
                                content = f.read()
                                f.seek(0)
                                f.truncate()
                                f.write(content.replace('doc/tex/makefile', ''))

                        # Run Configuration command
                        configcmd = "./configure" if package != 'metis' else ["make", "config"]
                        run_command(configcmd, currentdir=ext_dir)
                        try:
                            run_command("make", currentdir=ext_dir)
                        except:
                            raise Exception("{} installation has failed.".format(package))
                        # Save message to be shown at the end
                        if os.path.exists(ext_dir+ "/" + package):
                            final_message.append("{} has been successfully installed.".format(package))

                else:
                    # Clone from repo
                    git_clone(external_urls[package], branch)
                    # Save message to be shown at the end
                    final_message.append("{} has been sucessfully cloned.".format(package))
            else:
                if tarball:
                    final_message.append("{} has been already installed.".format(package))
                else:
                    # Checkout to the requested branch
                    os.chdir(top_dir + '/' + package)
                    subprocess.Popen(["git", "checkout", branch])
                    # Save message to be shown at the end
                    final_message.append("-- Skip cloning {}, because the folder already exists.".format(package))
                    final_message.append("-- Checking out {} ".format(package) + branch)
                    continue

        # Save post installation message if there is any.
        if package in messages.keys():
            final_message.extend(messages[package])

        # Change to top_dir
        os.chdir(top_dir)

    # Save post installation message about dunecontrol if need be.
    if not cleanup and any(x in pkg for pkg in packages for x in ["dumux","dune","opm"]):
        final_message.append("\n\nPlease run the following command (can be copied to command line):\n\n  ./dune-common/bin/dunecontrol --opts=./dumux/cmake.opts all")

    # If cleanup and only logfile in the external directory, remove the directory
    if os.path.isdir(ext_dir):
        _, _, files = next(os.walk(ext_dir))
        if cleanup and len(files)==1 and 'installexternal.log' in files:
            shutil.rmtree(ext_dir)

    return '\n'.join(final_message)

#################################################################
#################################################################
## (1/3) Define th necessary packages and their urls
#################################################################
#################################################################
dune_git_baseurl = "https://gitlab.dune-project.org/"
dumux_git_baseurl = "https://git.iws.uni-stuttgart.de/dumux-repositories/"
external_urls = {
    "dumux-lecture": dumux_git_baseurl + "dumux-lecture.git",
    "dumux-course": dumux_git_baseurl + "dumux-course.git",
    "dune-uggrid": dune_git_baseurl + "/staging/dune-uggrid.git",
    "dune-alugrid": dune_git_baseurl + "extensions/dune-alugrid.git",
    "dune-foamgrid": dune_git_baseurl + "extensions/dune-foamgrid.git",
    "dune-subgrid": "https://git.imp.fu-berlin.de/agnumpde/dune-subgrid.git",
    "dune-spgrid": dune_git_baseurl + "extensions/dune-spgrid.git",
    "dune-mmesh": dune_git_baseurl + "samuel.burbulla/dune-mmesh.git",
    "dune-functions": dune_git_baseurl + "staging/dune-functions.git",
    "dune-typetree": dune_git_baseurl + "staging/dune-typetree.git",
    "glpk": "http://ftp.gnu.org/gnu/glpk/glpk-4.60.tar.gz",
    "nlopt": "http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz",
    "opm-common": "https://github.com/OPM/opm-common",
    "opm-grid": "https://github.com/OPM/opm-grid",
    "metis": "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz",
    "gstat": "http://gstat.org/gstat.tar.gz",
}

packagenames = {
    "dumux-extensions": ["dumux-lecture", "dumux-course"],
    "dune-extensions": ["dune-uggrid", "dune-alugrid", "dune-foamgrid",
                        "dune-subgrid", "dune-spgrid", "dune-mmesh",
                        "dune-functions", "dune-typetree"],
    "functions": ["dune-functions", "dune-typetree"],
    "optimization": ["glpk", "nlopt"],
    "others": ["opm-common", "opm-grid", "metis", "gstat"]
}

messages ={
    'glpk': ["In addition, it might be necessary to set manually",
            "the glpk path in the CMAKE_FLAGS section of the .opts-file:",
            "  -DGLPK_ROOT=/path/to/glpk \\"],
    'dune-mmesh': ["Maybe you also have to install CGAL",
            "(see cgal.org/download.html)", "Maybe you also need to change your core dune modules' branches to their newest versions."],
    'opm-common': ["In addition, it might be necessary to set manually some",
                "CMake variables in the CMAKE_FLAGS section of the .opts-file:",
                "  -DUSE_MPI=ON",
                "Currently, compiling opm with clang is not possible.", "",
                "Maybe you also have to install the following packages (see the",
                " opm prerequisites at opm-project.org):",
                "  BLAS, LAPACK, Boost, SuiteSparse, Zoltan"]
}


#################################################################
#################################################################
## (2/3) Download/Config/Clean the requested packages
#################################################################
#################################################################
# Start download/configuration/cleaning tasks
final_message = install_external(args)

#################################################################
#################################################################
## (3/3) Show the final message
#################################################################
#################################################################
# Show final message
show_message(final_message)
