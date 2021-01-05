#!/usr/bin/env python3
import os
import sys
import string
import shutil
import argparse
import subprocess
from getmoduleinfo import getModuleFile
from getmoduleinfo import extractModuleInfos

# require python 3
if sys.version_info[0] < 3:
    sys.exit("\nERROR: Python3 required")

# input argument parser
parser = argparse.ArgumentParser(description="This script creates a docker image for a given module and module installation script.")
parser.add_argument('-m', '--modulepath', required=True, help='the path to the your module')
parser.add_argument('-i', '--installScript', required=True, help="Specify the installation script")
parser.add_argument('-t', '--templateFolder', required=False, help="Specify the folder with the template files")
args = vars(parser.parse_args());

# get information on the module
modulePath = args['modulepath']
moduleFolder = os.path.relpath(modulePath, os.path.join(modulePath, '../'))

modInfo = extractModuleInfos(getModuleFile(moduleFolder), ['Module', 'Maintainer'])
moduleName = modInfo['Module']
moduleMaintainer = modInfo['Maintainer']
dockerTag = moduleName.lower() # docker only supports lower case

# get folder with the template files
templateFolder = args['templateFolder']
if not templateFolder:
    templateFolder = os.path.join(modulePath, '../dumux/docker')
if not os.path.exists(templateFolder):
    sys.exit("Template file folder {} could not be found".format(templateFolder))

print("*"*54)
print("\n-- Creating a Docker image for module " + moduleName + " --\n" + "*"*54)

if os.path.exists("docker"):
    print("\nA docker folder already exists. Continue anyway? - will be overwritten - [y/N]\n")
    delete = input()
    if delete == "y" or delete == "Y":
        shutil.rmtree("docker")
        print("--> Deleted old docker folder.")
    else:
        sys.exit("Abort.")

os.mkdir("docker")
print("--> Created the folder 'docker'.")

# copy install script into docker folder and make it executable
installScriptArg = args['installScript']
installScriptName = os.path.split(installScriptArg)[1]
installScript = os.path.join(os.path.join(os.getcwd(), 'docker'), installScriptName)
shutil.copy(installScriptArg, installScript)
os.system("chmod +x {}".format(installScript))
print("--> Using install script: {} to install dependencies for module {}.".format(installScript, moduleName))

# substitute content from template and write to target
def substituteAndWrite(template, target, mapping):
    if not os.path.exists(template):
        sys.exit("Template file '" + template + "' could not be found")
    with open(target, 'w') as targetFile:
        raw = string.Template(open(template).read())
        targetFile.write(raw.substitute(**mapping))

# write setpermissions helper script
template = os.path.join(templateFolder, 'setpermissions.sh.template')
target = os.path.join(os.getcwd(), 'docker/setpermissions.sh')
substituteAndWrite(template, target, {})
print("--> Created permission helper script for easier container setup.")

# write welcome message file
template = os.path.join(templateFolder, 'WELCOME.template')
target = os.path.join(os.getcwd(), 'docker/WELCOME')
substituteAndWrite(template, target, {'modName': moduleName, 'modFolder': moduleFolder})
print("--> Created welcome message displayed on Docker container startup.")

# write readme file
template = os.path.join(templateFolder, 'README.md.template')
target = os.path.join(os.getcwd(), 'docker/README.md')
substituteAndWrite(template, target, {'modName': moduleName, 'dockerTag': dockerTag})
print("--> Created README.md on how to use the docker image.")

# write helper file for container spin-up (make it executable after creation)
template = os.path.join(templateFolder, 'docker.sh.template')
target = os.path.join(os.getcwd(), 'docker/docker_{}.sh'.format(dockerTag))
substituteAndWrite(template, target, {'dockerTag': dockerTag})
os.system("chmod +x " + target)
print("--> Created helper script to spin up the docker container.")

# write the docker file
template = os.path.join(templateFolder, 'Dockerfile.template')
target = os.path.join(os.getcwd(), 'docker/Dockerfile')
substituteAndWrite(template, target,
                  {
                      'modName': moduleName,
                      'modMaintainer': moduleMaintainer,
                      'dockerTag': dockerTag,
                      'instScript': installScriptName
                  })
print("--> Created Dockerfile. You can adapt it to your needs.")
print()
print("Do you want to directly build the Docker image? [y/N]")

build = input()
if build == "y" or build == "Y":
    print("Building Docker image... this may take several minutes.")
    try:
        os.chdir('docker')
        subprocess.run(['docker', 'build', '-f', 'Dockerfile', '-t', dockerTag, '.'], check=True)
        os.chdir('../')
    except:
        os.chdir('../')
        sys.exit("ERROR: docker image build failed")
    print()
    print("Successfully built docker image: " + dockerTag + ". Have a look at docker/README.md.")
    print("Check the container running docker run -it " + dockerTag + " /bin/bash in the same")
    print("directory as the Dockerfile and try using the convenience script docker_" + dockerTag + ".sh")
    print("See docker/README.md for more information.")
else:
    print("You can build your Docker image later by running docker build -f Dockerfile -t " + dockerTag)
    print("from within the folder 'docker' that was created by this script, and where the Dockerfile is.")
