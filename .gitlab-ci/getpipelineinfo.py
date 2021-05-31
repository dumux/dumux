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


def performApiQuery(command, err='API query unsuccessful'):
    return runCommand(command, suppressTraceBack=True, errorMessage=err)


def getPipeLinesApiURL(apiURL):
    return apiURL.rstrip('/') + '/pipelines/'


def getPipelines(apiURL, token, filter=''):
    queryURL = getPipeLinesApiURL(apiURL) + filter
    queryCmd = 'curl --header "token={}" "{}"'.format(token, queryURL)
    pl = performApiQuery(queryCmd, 'Could not retrieve pipelines')
    return json.loads(pl)


def getPipelineInfo(apiURL, token, pipeLineId, infoString):
    queryURL = getPipeLinesApiURL(apiURL) + str(pipeLineId) + '/' + infoString
    queryCmd = 'curl --header "token={}" "{}"'.format(token, queryURL)
    pl = performApiQuery(queryCmd, 'Could not retrieve pipeline info')
    return json.loads(pl)


def getPipeLineJobs(apiURL, token, pipeLineId):
    return getPipelineInfo(apiURL, token, pipeLineId, 'jobs')


def findPipeline(pipeLines, predicate):
    for pipeLine in pipeLines:
        if predicate(pipeLine):
            return pipeLine
    return None


parser = ArgumentParser(
    description='Find and print information on a previously run pipeline'
)

parser.add_argument('-p', '--print-format',
                    required=True, choices=['pipeline-id', 'commit-sha'],
                    help='Switch between reporting the pipeline-id/commit-sha')
parser.add_argument('-l', '--look-for',
                    required=True, choices=['latest', 'HEAD'],
                    help='Define how to search for pipelines')
parser.add_argument('-t', '--access-token',
                    required=True,
                    help='The token to post read requests to the GitLab API')
parser.add_argument('-u', '--project-api-url',
                    required=False,
                    default='https://git.iws.uni-stuttgart.de/api/v4/projects/31/',
                    help='The token to post read requests to the GitLab API')
parser.add_argument('-f', '--filter',
                    required=False, default='?status=success',
                    help='Pipeline query filter (default: "?status=success"')
parser.add_argument('-e', '--exclude-jobs',
                    required=False, nargs='*',
                    default=['check-pipeline-status'],
                    help='Exclude pipelines that contain the given jobs')
args = vars(parser.parse_args())

apiURL = args['project_api_url']
token = args['access_token']
pipeLines = getPipelines(apiURL, token, args['filter'])

currentBranch = runCommand('git branch --show-current').strip('\n')
currentCommitInfo = runCommand('git show HEAD').split('\n')
headIsMergeCommit = 'Merge:' in currentCommitInfo[1]

headSHA = runCommand('git rev-list HEAD --max-count=1').strip('\n')
if headIsMergeCommit:
    preSHA = runCommand('git rev-list HEAD --max-count=2').split('\n')[1]
    mrSHA = currentCommitInfo[1].split()[2]


def checkBranch(pipeLine):
    return pipeLine['ref'] == currentBranch


def checkCommit(pipeLine):
    sha = pipeLine['sha']
    if headIsMergeCommit:
        return sha == headSHA or sha == preSHA or mrSHA in sha
    return sha == headSHA


def skip(pipeLine):
    jobs = getPipeLineJobs(apiURL, token, pipeLine['id'])
    jobNames = [j['name'] for j in jobs]
    return any(j in jobNames for j in args['exclude_jobs'])


if args['look_for'] == 'HEAD':
    pipeLine = findPipeline(
        pipeLines, lambda p: checkCommit(p) and not skip(p)
    )
elif args['look_for'] == 'latest':
    pipeLine = findPipeline(
        pipeLines, lambda p: checkBranch(p) and not skip(p)
    )

if pipeLine is not None:
    if args['print_format'] == 'pipeline-id':
        print(pipeLine['id'])
    elif args['print_format'] == 'commit-sha':
        print(pipeLine['sha'])
else:
    sys.exit('Could not find a succesful pipeline')
