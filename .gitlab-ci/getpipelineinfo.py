import json
import sys
import os
from argparse import ArgumentParser

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))
    from common import runCommand
except Exception:
    sys.exit('Could not import common module')


def findPipeline(pipelineFile, predicate):
    with open(pipelineFile) as pipeLines:
        for pipeLine in json.load(pipeLines):
            if predicate(pipeLine):
                return pipeLine
    return None


parser = ArgumentParser(
    description='Find and print information on a previously run pipeline'
)

parser.add_argument('-s', '--status',
                    required=False, default='success',
                    help='Define the status of the pipeline to be found')
parser.add_argument('-p', '--print-format',
                    required=True, choices=['pipeline-id', 'commit-sha'],
                    help='Switch between reporting the pipeline-id/commit-sha')
parser.add_argument('-l', '--look-for',
                    required=True, choices=['latest', 'HEAD'],
                    help='Define how to search for pipelines')
parser.add_argument('-f', '--status-file',
                    required=True,
                    help='The .json file with the latest pipeline statuses')
args = vars(parser.parse_args())

currentBranch = runCommand('git branch --show-current').strip('\n')
currentCommitInfo = runCommand('git show HEAD').split('\n')
headIsMergeCommit = 'Merge:' in currentCommitInfo[1]

headSHA = runCommand('git rev-list HEAD --max-count=1').strip('\n')
if headIsMergeCommit:
    preSHA = runCommand('git rev-list HEAD --max-count=2').split('\n')[1]
    mrSha = currentCommitInfo[1].split()[2]


def checkCommitSHA(pipeLine):
    sha = pipeLine['sha']
    if headIsMergeCommit:
        return sha == headSHA or sha == preSHA or mrSha in sha
    return sha == headSHA


def checkStatus(pipeLine): return pipeLine['status'] == args['status']
def checkBranch(pipeLine): return pipeLine['ref'] == currentBranch
def find(predicate): return findPipeline(args['status_file'], predicate)


if args['look_for'] == 'HEAD':
    pipeLine = find(lambda p: checkStatus(p) and checkCommitSHA(p))
elif args['look_for'] == 'latest':
    pipeLine = find(lambda p: checkStatus(p) and checkBranch(p))

if pipeLine is not None:
    if args['print_format'] == 'pipeline-id':
        print(pipeLine['id'])
    elif args['print_format'] == 'commit-sha':
        print(pipeLine['sha'])
else:
    sys.exit('Could not find a succesful pipeline')
