"""Local-gap input validation, legacy defaults, and Python/C++ parameter wiring."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('run_graphite_local_adhesion',
                                            ROOT / 'slurry/drivers/gr_re2/run_graphite.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def case():
    return json.loads((ROOT / 'slurry/cases/pure_gr.json').read_text())


def invalid_cases():
    cases = []
    for key, value in [('local_gap_fraction', -0.1), ('local_gap_fraction', 1.01),
                       ('local_gap_fraction', math.nan), ('local_gap_fraction', math.inf),
                       ('local_gap_m', 0), ('local_gap_m', 3e-9),
                       ('local_switch_excess_gap_m', -1e-9),
                       ('local_switch_excess_gap_m', 1e-8),
                       ('local_cutoff_excess_gap_m', 0),
                       ('local_cutoff_excess_gap_m', 4e-7)]:
        config = case()
        config['interaction'][key] = value
        cases.append((key + '=' + str(value), config))
    config = case()
    config['rough_contact']['enabled'] = False
    cases.append(('rough_contact disabled', config))
    return cases


class LocalAdhesionInputTests(unittest.TestCase):
    def test_case_round_trip_and_recorded_parameters(self):
        cfg, metadata = RUNNER.resolve(case())
        values = RUNNER.solver_values(cfg, Path('/tmp/run'), Path('/tmp/particles.csv'), 0)
        self.assertEqual(values['sigma_lj'], 4.197e-10)
        self.assertEqual(values['local_gap'], 3e-10)
        self.assertEqual(values['local_gap_fraction'], 1)
        self.assertEqual(values['local_switch_excess_gap'], 2e-9)
        self.assertEqual(values['local_cutoff_excess_gap'], 1e-8)
        self.assertEqual(values['roughness_gap'], 2e-9)
        self.assertEqual(values['tangential_stiffness'], 80)
        self.assertEqual(values['sliding_friction'], 1)
        self.assertEqual(metadata['interaction'], cfg['interaction'])
        self.assertEqual(metadata['rough_contact'], cfg['rough_contact'])

    def test_missing_local_inputs_preserve_old_physics(self):
        original = case()
        original['interaction'] = {key: value for key, value in original['interaction'].items()
                                   if not key.startswith('local_')}
        original['interaction']['sigma_lj_m'] = 3e-9
        original['rough_contact']['tangential_stiffness_N_m'] = 9
        snapshot = copy.deepcopy(original)
        cfg, _ = RUNNER.resolve(original)
        self.assertEqual(cfg['interaction']['local_gap_fraction'], 0)
        self.assertEqual(cfg['interaction']['sigma_lj_m'], 3e-9)
        self.assertEqual(cfg['rough_contact']['tangential_stiffness_N_m'], 9)
        self.assertEqual(original, snapshot)

    def test_invalid_local_inputs(self):
        for name, cfg in invalid_cases():
            with self.subTest(name=name), self.assertRaises(ValueError):
                RUNNER.resolve(cfg)

    def test_disabled_local_correction_does_not_require_rough_contact(self):
        cfg = case()
        cfg['rough_contact']['enabled'] = False
        cfg['interaction']['local_gap_fraction'] = 0
        RUNNER.resolve(cfg)

    def test_old_binary_only_rejected_for_active_local_correction(self):
        cfg, _ = RUNNER.resolve(case())
        with self.assertRaisesRegex(ValueError, 'Rebuild'):
            RUNNER.require_local_adhesion_build(cfg, {'rough_contact': True})
        RUNNER.require_local_adhesion_build(cfg, {'local_gap_adhesion': True})
        cfg['interaction']['local_gap_fraction'] = 0
        RUNNER.require_local_adhesion_build(cfg, {})


@unittest.skipUnless(shutil.which('g++'), 'C++ parser integration requires g++')
class CppLocalAdhesionInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='local-adhesion-config-')
        cls.directory = Path(cls.temporary.name)
        source = cls.directory / 'parse.cpp'
        source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  try {
    const auto c=slurry::gr_re2::parseConfig(argc,argv);
    std::cout<<std::setprecision(17)<<c.sigma_lj<<" "<<c.local_gap<<" "
      <<c.local_gap_fraction<<" "<<c.local_switch_excess_gap<<" "
      <<c.local_cutoff_excess_gap<<" "<<c.roughness_gap<<" "<<c.tangential_stiffness;
  } catch(const std::exception& e) { std::cerr<<e.what();return 2; }
}
''')
        cls.executable = cls.directory / 'parse'
        subprocess.run(['g++', '-std=c++17', '-I' + str(ROOT / 'olb-1.9r0/src/slurry/gr_re2'),
                        str(source), '-o', str(cls.executable)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def parse(self, cfg, remove_local=False):
        # The invalid cases deliberately bypass Python resolve to exercise C++.
        baseline, _ = RUNNER.resolve(case())
        baseline['interaction'].update(cfg['interaction'])
        baseline['rough_contact'].update(cfg['rough_contact'])
        values = RUNNER.solver_values(baseline, self.directory, self.directory / 'particles.csv', 0)
        if remove_local:
            values = {key: value for key, value in values.items() if not key.startswith('local_')}
        path = self.directory / 'run.cfg'
        path.write_text(''.join('{}={}\n'.format(key, value) for key, value in values.items()))
        return subprocess.run([str(self.executable), '--config', str(path)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    def test_generated_parameters_reach_cpp_unchanged(self):
        result = self.parse(case())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([float(value) for value in result.stdout.split()],
                         [4.197e-10, 3e-10, 1, 2e-9, 1e-8, 2e-9, 80])

    def test_cpp_rejects_invalid_local_inputs(self):
        for name, cfg in invalid_cases():
            with self.subTest(name=name):
                result = self.parse(cfg)
                self.assertEqual(result.returncode, 2, result.stdout)

    def test_cpp_missing_fraction_remains_disabled(self):
        result = self.parse(case(), remove_local=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(float(result.stdout.split()[2]), 0)


if __name__ == '__main__':
    unittest.main()
