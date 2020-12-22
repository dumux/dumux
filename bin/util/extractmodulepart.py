#!/usr/bin/env python3
import sys, glob, os, subprocess
import argparse
import shutil
import os.path

"""
This is a python script for extracting modules.
It was originally written to solve the problem that the bash shell script could only run in a specific environment.
So the whole design philosophy is to maximize the portability, i.e. write in a cross-platform style.
"""

# def somes function to print the help message, explanations, guidience.

def show_helpmesseage():
    print("\n""USAGE: "+os.path.basename(__file__)+" module_dir FOLDER_1 [FOLDER_2 ...]\n\n"
          "module_dir is the folder containing the DUNE module from which you\n"
          "want to extract. The script has to be called one level above it.\n\n"
          "The FOLDERs need to indicate subfolders of module_dir. At least one\n"
          "of them has to contain a source file *.cc of an executable for which\n"
          "you would like to timber a table in dumux-pub.")
    exit(1)

# check if help is needed
if (len(sys.argv) < 2):
    show_helpmesseage()
if (str(sys.argv[1]) == "--help" or
    str(sys.argv[1]) == "-help" or
    str(sys.argv[1]) == "help"):
    show_helpmesseage()
module_dir = str(sys.argv[1])

# if module_dir contains a slash as last character, delete it
if module_dir.endswith('/'):
   module_dir = module_dir[:-1]

# check if we are above module_dir
if not (os.path.isdir(module_dir)):
    print("ERROR: you need to run the script one level above the folder "+module_dir+".")
    print("Run \""+os.path.basename(__file__)+" --help\" for details.")
    exit(1)

# determine all source files in the paths passed as arguments
os.chdir(module_dir)
module_full_path=os.getcwd()
all_sources=[]
all_directories=[]
for dir_path in sys.argv[2:]:
    stripped_path = dir_path.removeprefix(module_dir)
    directories = " " + stripped_path
    all_directories.append(stripped_path)
    os.chdir( os.path.join(module_full_path + stripped_path))
    for file in glob.glob("*.cc"):
        sources = os.path.join(stripped_path, file)
        all_sources.append(sources)
os.chdir(module_full_path)
os.chdir("..")

# check if sources have been obtained
for source in all_sources:
    contracted = str(os.popen("echo \"" + source + "\" | tr -d \" tnr\"").read().rstrip("\n"))
    if ( contracted == "" ):
        print("ERROR: no source files *.cc found in the directories "+sys.argv[2]+".")
        print("Be sure to provide a list of paths as arguments to this script.")
        print("Run \""+os.path.basename(__file__)+" --help\" for details.")
        exit(1)

# try to find the duneproject script
dune_project = shutil.which('duneproject', path = "dune-common/bin")
if (dune_project == None):
    print("ERROR: Could not find duneproject.")
    print("Be sure to either have duneproject in your search path")
    print("or to run this script from a directory that contains duneproject.")
    exit(1)
else:
    print(dune_project)

# give explanations
print("\n""This script will\n"
      "- extract the following sub-folders of "+ module_dir +":\n\n")
for dir_path in all_directories:
    print("  "+ dir_path + ",")
print("\n""  and all headers in "+ module_dir + " that are required to build the\n"
      "  executables from the sources\n\n")
for source in all_sources:
    print("  "+ source +",")
print("\n- copy the extracted files into a freshly created DUNE module, retaining the\n"
      "  directory structure,\n\n"
      "- update/create all required files like CMakeLists.txt,\n\n"
      "- store the versions of all used Dune module\n\n"
      "- and extract their modifications as patches.\n\n"
      "Thus, you receive a fully-working DUNE module containing the subset of\n"
       + module_dir + " that is required to run your application.\n"
      "duneproject will be run now. The new module should NOT depend on the\n"
      "module in "+ module_dir + ".\n\n"
      "Read the above and press [Enter] to proceed...",end="")

# run duneproject
old_ls = os.listdir()
subprocess.call([dune_project], shell=True)
new_ls  = os.listdir()

# determine the new module/directory name
module_name = (set(new_ls) - set(old_ls)).pop()
if ( module_name  == "" ):
    print("ERROR: could not find new module. Aborting.")
    exit(1)
else:
    print()
    print( os.path.basename(__file__) + ": Found new module " + module_name)
print("Determining required headers...",end="")
os.chdir(module_name)

# extract all headers
import re
import threading

all_file = []
def search_headers(c_file):
    global all_file
    f = open(c_file,'r')
    content = f.read()
    f.close()
    all_info = re.findall(r'(?<=#include <).+?(?=>)',content)
    # all_info_problem = re.findall(r'(?<=#include ").+?(?=")',content)
    # all_info = all_info_dumux + all_info_problem
    if all_info == []:
        sys.exit(1)
    else:
        for one_info in all_info:
            if one_info.startswith("dumux"):
                b = os.path.join(module_full_path, one_info)
            else:
                continue
            all_file.append(b)
            thread_ = threading.Thread(target = search_headers, args = (b, ))
            thread_.start()

for source in all_sources:
    source_abspath = module_full_path+source
    search_headers(source_abspath)
    print(all_file)
    for i in all_file:
        dir_path = os.path.dirname(os.path.realpath(i)).replace("dumux",module_name,1)
        os.makedirs(dir_path, exist_ok=True)
        shutil.copy(i,dir_path)
    source_dir = source_abspath.removesuffix("main.cc")
    source_path = source_dir.replace("dumux",module_name,1)
    # os.makedirs(source_path, exist_ok=True)
    shutil.copytree(source_dir, source_path)

# delete all architecture-dependent files and unneeded directories
rmfilenames = ["Makefile.in", "Makefile", '*.o', '*.deps/*']
for rmfile in rmfilenames:
    for filename in glob.glob(rmfile):
        os.remove(filename)

shutil.rmtree('dune')
shutil.rmtree('src')

# set CMake File for each directory

# move patches folder into module if existing
if (os.path.isdir("patches") ):
    subprocess.call(["mv","patches",module_name],shell=True)

# output guidence for users
print("\n"+"*"*80+"\n"
      "The extracted module is contained in the subfolder \""+module_name+"\".\n"
      "You can build it using \"dunecontrol ... --only="+ module_name +" all\".\n"
      +"*"*80+"\n"
      "BEFORE building, you can add the module to dumux-pub by something like:\n"
      "(Rename module name if it does not match the AuthorLastNameYearx scheme\n"
      "and commit it to the Git repository dumux-pub)\n"
      "git clone https://git.iws.uni-stuttgart.de/dumux-pub/AuthorLastNameYearx.git \n"
      "sed -i '/Module:/c\Module: AuthorLastNameYearx' "+ module_name+"/dune.module\n"
      "mv "+ module_name +"/* dumux-pub/AuthorLastNameYearx/.\n"
      "cd AuthorLastNameYearx\n"
      "git commit -a\n"
      "git push")
