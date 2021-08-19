#!/usr/bin/env python3

from statistics import mean
from math import log, isnan, isinf
from argparse import ArgumentParser
import subprocess


class Test:
    def __init__(self, args):
        self.args = args
        self.darcySuffix = 'darcy'
        self.freeFlowSuffix = 'freeFlow'
        self.numRefinements = 2

    def run(self):
        self._removeLogFile(self.darcySuffix)
        self._removeLogFile(self.freeFlowSuffix)

        self._compile()
        for i in range(0, self.numRefinements+1):
            self._runDumux(i)

        darcyLog = self._getLogFileName(self.darcySuffix)
        freeFlowLog = self._getLogFileName(self.freeFlowSuffix)

        errorsDarcy = self._collectErrors(darcyLog)
        ratesDarcy = self._collectRates(errorsDarcy)
        errorsFreeFlow = self._collectErrors(freeFlowLog)
        ratesFreeFlow = self._collectRates(errorsFreeFlow)

        self._checkMeanRates(ratesDarcy, ratesFreeFlow)

    def _checkMeanRates(self, ratesDarcy, ratesFreeFlow):
        def checkRate(rates, key, expected):
            m = mean(rates[key])
            if m < expected:
                raise Exception('Rate {} is below expected value {}: {}'
                                .format(key, expected, m))
            else:
                print(
                    'Test for "{}" successful: '
                    'computed vs expected rate {:.3f}/{:.3f}'
                    .format(key, m, expected)
                )

        expectedP = float(self.args['expectedpressurerate'])
        expectedV = float(self.args['expectedvelocityrate'])
        checkRate(ratesDarcy, 'p_rate', expectedP)
        checkRate(ratesFreeFlow, 'p_rate', expectedP)
        checkRate(ratesFreeFlow, 'vx_rate', expectedV)
        checkRate(ratesFreeFlow, 'vy_rate', expectedV)

    def _compile(self):
        subprocess.run(['make', self.args['exe']], check=True)

    def _runDumux(self, refIdx):
        command = [
            './{}'.format(self.args['exe']), self.args['inputfile'],
            '-Problem.TestCase', self.args['case'],
            '-Grid.Refinement', str(refIdx)
        ]
        if args['dumux_args']:
            command.extend(args['dumux_args'])
        subprocess.run(command, check=True)

    def _getLogFileName(self, suffix):
        return self.args['case'] + '_' + suffix + '.log'

    def _removeLogFile(self, suffix):
        logFileName = self._getLogFileName(suffix)
        subprocess.run(['rm', '-f', logFileName], check=True)
        print("Removed old log file ({})!".format(logFileName))

    def _collectErrors(self, logFileName):
        result = {}
        with open(logFileName, 'r') as logFile:
            for line in logFile:
                line = line.strip(r'[ConvergenceTest]')
                errors = line.split('L2(')
                errors = filter(lambda e: '=' in e, errors)
                errors = [err.split('=') for err in errors]
                for name, val in errors:
                    name = name.strip().strip(')')
                    val = val.strip()
                    if name not in result:
                        result[name] = []
                    result[name].append(float(val))
        return result

    def _collectRates(self, errors):
        rates = {}
        for name, errValues in errors.items():
            for i in range(len(errValues)-1):
                if isnan(errValues[i]) or isinf(errValues[i]):
                    continue
                if not (
                    (errValues[i] < 1e-12 or errValues[i+1] < 1e-12) and
                    (errValues[i] < 1e-12 or errValues[i+1] < 1e-12) and
                    (errValues[i] < 1e-12 or errValues[i+1] < 1e-12)
                ):
                    rate = (log(errValues[i]) - log(errValues[i+1]))/log(2)
                    key = name + '_rate'
                    if key not in rates:
                        rates[key] = []
                    rates[key].append(rate)
                else:
                    raise Exception("Error: exact solution!?")
        return rates


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-e', '--exe', required=True)
    parser.add_argument('-c', '--case', required=True)
    parser.add_argument('-rp', '--expectedpressurerate', default='1.75')
    parser.add_argument('-rv', '--expectedvelocityrate', default='1.75')
    parser.add_argument('-i', '--inputfile', default='params.input')

    args, unknown = parser.parse_known_args()
    args = vars(args)
    args['dumux_args'] = unknown

    Test(args).run()
