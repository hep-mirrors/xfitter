#!/usr/bin/env python3
"""Check the CERES EoE fit numerically; write a GitLab-compatible JUnit report."""
import argparse
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    result = float(value.replace('D', 'E'))
    require(math.isfinite(result), 'Non-finite number: ' + value)
    return result


class Validation:
    def __init__(self, run):
        self.run = run
        self.output = run / 'output'
        self.expected = json.loads(Path(__file__).with_name('expected.json').read_text())

    def results(self):
        return (self.output / 'Results.txt').read_text()

    def match(self, pattern, text):
        match = re.search(pattern, text, re.M)
        require(match is not None, 'Missing output: ' + pattern)
        return match

    def factors(self):
        log = (self.run / 'xfitter.log').read_text()
        gof = number(self.match(r'^Chi2 Bartlett factor.*=\s*(\S+)', log)[1])
        ci = number(self.match(r'^CI Bartlett factor.*=\s*(\S+)', log)[1])
        return gof, ci

    def sources(self):
        rows = {}
        for line in self.results().splitlines():
            fields = line.split()
            if len(fields) == 12 and fields[0].isdigit() and '@eps=' in fields[-1]:
                require(fields[1] not in rows, 'Duplicate EoE source: ' + fields[1])
                rows[fields[1]] = ([number(x) for x in fields[2:11]], fields[-1])
        return rows

    def convergence(self):
        require((self.output / 'Status.out').read_text().strip() == 'OK',
                'CERES did not report successful convergence')
        log = (self.run / 'xfitter.log').read_text()
        self.match(r'Termination:\s+CONVERGENCE', log)
        require('CERES minimisation' in log, 'CERES was not run')
        require('UNAVAILABLE' not in self.results(), 'Bartlett corrections unavailable')
        require(not re.search(r'W:.*(?:[Nn]ot converg|[Nn]on.converg|maximum.*iteration)', log),
                'Profiling reported a convergence warning')

    def chi2_and_bartlett(self):
        expected = self.expected
        text = self.results()
        ndf = expected['npoints'] - expected['npoi']
        raw = self.match(r'^Chi2 after minimisation\s+(\S+)\s+(\d+)\s+(\S+)', text)
        corrected = self.match(r'^Corrected Chi2\s+(\S+)\s+(\d+)\s+(\S+)', text)
        for row, key in [(raw, 'chi2'), (corrected, 'corrected_chi2')]:
            value = number(row[1])
            require(int(row[2]) == ndf, 'Incorrect final degrees of freedom')
            require(abs(number(row[3]) - value / ndf) < 0.0006, 'Incorrect chi2/ndf')
            require(abs(value - expected[key]) <= expected[key + '_tolerance'],
                    '{} outside regression tolerance: {}'.format(key, value))
        counts = self.match(r'with ndf_bart = npoints - nPOI =\s*(\d+)\s*,\s*nPOI =\s*(\d+)', text)
        require((int(counts[1]), int(counts[2])) == (ndf, expected['npoi']),
                'Bartlett parameter count is incorrect')
        gof, ci = self.factors()
        for value, key in [(gof, 'gof_factor'), (ci, 'ci_factor')]:
            require(abs(value - expected[key]) <= expected[key + '_tolerance'],
                    '{} outside regression tolerance: {}'.format(key, value))
        require(0 < gof < 1 and ci > 1, 'Missing or invalid Bartlett scaling')
        require(abs(number(corrected[1]) - number(raw[1]) * gof) < 0.015,
                'Corrected chi2 does not use the reported factor')
        rows = self.sources().values()
        require(abs(gof - 1 / (1 + sum(v[7] for v, _ in rows) / ndf)) < 2e-6,
                'GoF factor disagrees with source contributions')
        require(abs(ci - math.sqrt(1 + sum(v[8] for v, _ in self.sources().values()) /
                                  expected['npoi'])) < 2e-6,
                'CI factor disagrees with source contributions')

    def eoe_sources(self):
        expected = self.expected
        rows = self.sources()
        names = set(expected['external'] + expected['profiled'])
        require(set(rows) == names, 'Incorrect set of EoE sources')
        for name, (values, kind) in rows.items():
            prefix = ':E:' if name in expected['external'] else ':N:'
            require(kind.startswith(prefix), 'Wrong nuisance treatment: ' + name)
            require(abs(values[3] - expected['epsilon']) < 1e-8 and
                    abs(number(kind.split('@eps=')[1]) - expected['epsilon']) < 1e-8,
                    'Incorrect epsilon: ' + name)
            require(values[1] > 0 and values[2] > 0, 'Invalid nuisance errors: ' + name)
            require(abs(values[2] - values[1] * math.sqrt(1 + values[5])) < 0.0002,
                    'Incorrect corrected nuisance uncertainty: ' + name)

    def parameter_output(self):
        raw = {}
        for match in re.finditer(r'^  (\w+): \[ ([^,]+), ([^ ]+) \]',
                                 (self.output / 'pars.yaml').read_text(), re.M):
            raw[match[1]] = (number(match[2]), number(match[3]))
        rows = {}
        indices = []
        for line in (self.output / 'parsout_0').read_text().splitlines():
            match = self.match(r"^\s*(\d+)\s+'([^']+)'\s+(\S+)\s+(\S+)\s*$", line)
            indices.append(int(match[1]))
            require(match[2] not in rows, 'Duplicate fitted parameter')
            rows[match[2]] = (number(match[3]), number(match[4]))
        count = self.expected['npoi'] + len(self.expected['external'])
        require(len(rows) == count and set(rows) == set(raw), 'Missing fitted parameters')
        require(indices == list(range(count)), 'CERES parameter indices changed')
        sources = self.sources()
        _, ci = self.factors()
        for name, (value, error) in rows.items():
            require(abs(value - raw[name][0]) <= 2e-6, 'Fitted value not restored: ' + name)
            require(error > 0 and raw[name][1] > 0, 'Invalid parameter error: ' + name)
            if name in self.expected['external']:
                require(abs(error - sources[name][0][2]) < 0.00006,
                        'External error disagrees with Results.txt: ' + name)
            else:
                require(abs(error - raw[name][1] * ci) < 5e-6,
                        'PDF parameter error lacks CI scaling: ' + name)

    def covariance(self):
        text = (self.output / 'ceres.out.txt').read_text()
        lines = text.split('COVARIANCE MATRIX', 1)[1].strip().splitlines()
        names = lines[0].split()
        count = self.expected['npoi'] + len(self.expected['external'])
        require(len(names) == count, 'Wrong covariance dimension')
        for i, line in enumerate(lines[1:count + 1]):
            fields = line.split()
            require(len(fields) == count + 1 and fields[0] == names[i], 'Malformed covariance')
            values = [number(x) for x in fields[1:]]
            require(values[i] > 0, 'Non-positive parameter variance')
        require(len(lines) >= count + 1, 'Truncated covariance')

    def pdf_bands(self):
        count = self.expected['npoi'] + len(self.expected['external'])
        files = [self.output / 'pdfs_q2val_01.txt']
        bands = sorted(self.output.glob('pdfs_q2val_s*s_01.txt'))
        require({p.name for p in bands} ==
                {'pdfs_q2val_s{:02d}s_01.txt'.format(i) for i in range(1, count + 1)},
                'Missing or unexpected symmetric PDF bands')
        grid = None
        for path in files + bands:
            lines = path.read_text().splitlines()
            header = [number(x) for x in lines[0].split()]
            require(len(header) == 5 and abs(header[0] - 4) < 1e-8, 'Wrong PDF Q2')
            nrows, ncols = int(header[1]), int(header[2]) + 1
            require(nrows > 0 and len(lines) == nrows + 2, 'Truncated PDF table: ' + path.name)
            table = [[number(x) for x in line.split()] for line in lines[2:]]
            require(all(len(row) == ncols for row in table), 'Malformed PDF table: ' + path.name)
            xgrid = [row[0] for row in table]
            require(all(x < y for x, y in zip(xgrid, xgrid[1:])), 'Unordered PDF x grid')
            if grid is None:
                grid = xgrid
            require(xgrid == grid, 'Band x grid differs from central PDF')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    args = parser.parse_args()
    validation = Validation(args.run_directory)
    checks = ['convergence', 'chi2_and_bartlett', 'eoe_sources',
              'parameter_output', 'covariance', 'pdf_bands']
    suite = ET.Element('testsuite', name='errors-on-errors-CERES', tests=str(len(checks)))
    failures = 0
    for name in checks:
        case = ET.SubElement(suite, 'testcase', classname='CERES.EoE', name=name)
        try:
            getattr(validation, name)()
            print('PASS: ' + name)
        except (OSError, ValueError, IndexError, KeyError, ZeroDivisionError) as error:
            failures += 1
            ET.SubElement(case, 'failure', message=str(error)).text = str(error)
            print('FAIL: {}: {}'.format(name, error))
    suite.set('failures', str(failures))
    ET.ElementTree(suite).write(str(args.run_directory / 'validation.xml'), encoding='utf-8',
                               xml_declaration=True)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
