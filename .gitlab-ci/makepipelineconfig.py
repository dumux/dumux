#!/usr/bin/env python3

import os
import sys
import string
from argparse import ArgumentParser

# require Python 3
if sys.version_info.major < 3:
    sys.exit('Python 3 required')

parser = ArgumentParser(description='Generate dumux test pipeline .yml file')
parser.add_argument('-o', '--outfile',
                    required=True,
                    help='Specify the file to write the pipeline definition')
parser.add_argument('-a', '--affectedtestsonly',
                    required=False,
                    action='store_true',
                    help='Use this flag to create a pipeline that runs only '
                         'those tests that are affected by changes w.r.t to '
                         'the origin/master branch')
parser.add_argument('-t', '--template',
                    required=False,
                    default='.gitlab-ci/default.yml.template',
                    help='Specify the template .yml file to be used')
parser.add_argument('-i', '--indentation',
                    required=False,
                    default=4,
                    help='Specify the indentation for the script commands')
args = vars(parser.parse_args())


# substitute content from template and write to target
def substituteAndWrite(mapping):

    template = args['template']
    if not os.path.exists(template):
        sys.exit("Template file '" + template + "' could not be found")

    with open(args['outfile'], 'w') as ymlFile:
        raw = string.Template(open(template).read())
        ymlFile.write(raw.substitute(**mapping))


commandIndentation = ' '*args['indentation']
with open(args['outfile'], 'w') as ymlFile:

    def wrapDuneControl(command):
        return 'dunecontrol --opts=$DUNE_OPTS_FILE --current ' + command

    def makeYamlList(commands, indent=commandIndentation):
        commands = [indent + '- ' + comm for comm in commands]
        return '\n'.join(commands)

    # if no configuration is given, build and run all tests (skip select stage)
    if not args['affectedtestsonly']:
        buildCommand = [wrapDuneControl('all'),
                        wrapDuneControl('bexec make -k -j4 build_tests')]
        testCommand = [wrapDuneControl('bexec dune-ctest'
                                       ' -j4 --output-on-failure')]

        substituteAndWrite({'build_script': makeYamlList(buildCommand),
                            'test_script': makeYamlList(testCommand),
                            'stages': makeYamlList(['build', 'test'], '  '),
                            'test_select_job': '',
                            'build_needs': ''})

    # otherwise, add a stage that detects the tests to be run first
    else:
        selectStageName = 'configure'
        selectJobName = 'select tests'

        stages = makeYamlList([selectStageName, 'build', 'test'], '  ')
        buildNeeds = '\n'.join(['  needs:',
                                '    - job: {}'.format(selectJobName),
                                '      artifacts: true'])

        selectJob = '\n'.join([selectJobName + ':',
                               '  stage: {}'.format(selectStageName),
                               '  script:'])
        selectJob += '\n'
        selectJob += makeYamlList([wrapDuneControl('all'),
                                   'pushd build-cmake',
                                   'python3 ../bin/testing/findtests.py'
                                   ' -f ../affectedtests.json'
                                   ' -t origin/master',
                                   'popd'])
        selectJob += '\n'
        selectJob += '\n'.join(['  artifacts:',
                                '    paths:',
                                '      - affectedtests.json',
                                '    expire_in: 3 hours'])

        buildCommand = [wrapDuneControl('all'),
                        'cp affectedtests.json build-cmake',
                        'pushd build-cmake',
                        'python3 ../bin/testing/runselectedtests.py '
                        ' -c affectedtests.json -b',
                        'popd']
        testCommand = [wrapDuneControl('all'),
                       'pushd build-cmake',
                       'python3 ../bin/testing/runselectedtests.py '
                       ' -c affectedtests.json -t',
                       'popd']

        substituteAndWrite({'build_script': makeYamlList(buildCommand),
                            'test_script': makeYamlList(testCommand),
                            'stages': stages,
                            'test_select_job': selectJob,
                            'build_needs': buildNeeds})
