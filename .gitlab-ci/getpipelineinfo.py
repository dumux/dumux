# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

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


class APIRequester:
    def __init__(self, projectURL, token):
        self.projectURL = projectURL.rstrip('/')
        self.token = token

    def __call__(self, urlSuffix, msg='API request failed'):
        reqURL = self.projectURL + '/' + urlSuffix.lstrip('/')
        reqCmd = 'curl --header "token={}" "{}"'.format(self.token, reqURL)
        data = runCommand(reqCmd, suppressTraceBack=True, errorMessage=msg)
        return json.loads(data)


# if no id is given this returns all pipelines (// is fine)
def pipelineApiSuffix(id=''):
    return 'pipelines/' + str(id) + '/'


# if no commit sha is given this returns all commits (// is fine)
def commitApiSuffix(sha=''):
    return 'repository/commits/' + str(sha) + '/'


parser = ArgumentParser(
    description='Find and print information on a previously run pipeline'
)

parser.add_argument('-p', '--print-format',
                    required=True, choices=['pipeline-id', 'commit-sha'],
                    help='Switch between reporting the pipeline-id/commit-sha')
parser.add_argument('-l', '--look-for',
                    required=True, choices=['latest', 'latest-merge', 'HEAD'],
                    help='Define how to search for pipelines')
parser.add_argument('-m', '--max-tree-depth',
                    required=False, type=int, default=50,
                    help='Maximum number of revisions to consider'
                         ' in search for pipeline candidates')
parser.add_argument('-t', '--access-token',
                    required=True,
                    help='The token to post read requests to the GitLab API')
parser.add_argument('-u', '--project-api-url',
                    required=False,
                    default='https://git.iws.uni-stuttgart.de/api/v4/projects/31/',
                    help='The token to post read requests to the GitLab API')
parser.add_argument('-s', '--pipeline-status',
                    required=False, default='success',
                    help='Status of pipeline candidates (default: success')
parser.add_argument('-e', '--exclude-jobs',
                    required=False, nargs='*',
                    default=['check-pipeline-status'],
                    help='Exclude pipelines that contain the given jobs')
args = vars(parser.parse_args())


requester = APIRequester(args['project_api_url'], args['access_token'])


def isMergeCommit(commitSHA):
    commitInfo = runCommand(f'git show --no-patch {commitSHA}').split('\n')
    return commitInfo[1].startswith('Merge:')


def getLastPipeline(commitSHA):
    return requester(commitApiSuffix(commitSHA))['last_pipeline']


def hasMatchingStatus(pipeline):
    return pipeline['status'] == args['pipeline_status']


def hasExcludeJob(pipeLine):
    suffix = pipelineApiSuffix(pipeLine['id'])
    jobs = requester(suffix.rstrip('/') + '/jobs/')
    return any(j in jobs for j in args['exclude_jobs'])


def findPipeline(commits):
    for commit in commits:
        pipeLine = getLastPipeline(commit)
        if pipeLine is not None:
            if hasMatchingStatus(pipeLine) and not hasExcludeJob(pipeLine):
                return pipeLine
    return None


if args['look_for'] == 'HEAD':
    commits = [runCommand('git rev-list HEAD --max-count=1').strip('\n')]
    if isMergeCommit(commits[0]):
        preSHA = runCommand('git rev-list HEAD --max-count=2').split('\n')[1]
        commits.append(preSHA)
    pipeLine = findPipeline(commits)

elif args['look_for'] == 'latest':
    count = args['max_tree_depth']
    commits = filter(None, runCommand(f'git rev-list HEAD --max-count={count}').split('\n'))
    pipeLine = findPipeline(commits)

elif args['look_for'] == 'latest-merge':
    count = args['max_tree_depth']
    commits = filter(None, runCommand(f'git rev-list HEAD --merges --max-count={count}').split('\n'))
    pipeLine = findPipeline(commits)

if pipeLine is not None:
    if args['print_format'] == 'pipeline-id':
        print(pipeLine['id'])
    elif args['print_format'] == 'commit-sha':
        print(pipeLine['sha'])
else:
    sys.exit('Could not find a successful pipeline')
