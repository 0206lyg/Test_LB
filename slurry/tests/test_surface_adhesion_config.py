"""Surface adhesion opt-in, input boundaries, and Python-to-C++ physics wiring."""
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
SPEC = importlib.util.spec_from_file_location('surface_adhesion_runner', ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def case():
    return json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())


def invalid_cases():
    cases = []
    values = [('surface_adhesion', 1), ('surface_adhesion', 'true'),
              ('adhesion_work_J_m2', 0), ('adhesion_work_J_m2', 1e-5),
              ('adhesion_work_J_m2', math.nan), ('adhesion_work_J_m2', math.inf),
              ('adhesion_range_m', 0), ('adhesion_range_m', math.inf),
              ('adhesion_range_m', 4e-9), ('curvature_switch_gap_m', 2e-9),
              ('curvature_switch_gap_m', 2e-8), ('curvature_cutoff_gap_m', 5e-7),
              ('curvature_cutoff_gap_m', math.nan), ('sigma_lj_m', 4e-9)]
    for key, value in values:
        cfg = case()
        cfg['interaction'][key] = value
        cases.append((key+'='+str(value), cfg))
    for key, value in RUNNER.LOCAL_ADHESION_DEFAULTS.items():
        cfg = case()
        cfg['interaction'][key] = value
        cases.append(('mixed '+key, cfg))
    cfg = case()
    cfg['rough_contact']['enabled'] = False
    cases.append(('rough contact disabled', cfg))
    cfg = case()
    cfg['interaction']['surface_adhesion'] = False
    cases.append(('inactive surface inputs', cfg))
    return cases


class SurfaceAdhesionInputTests(unittest.TestCase):
    def test_resolve_preserves_input_and_all_physical_fields(self):
        original = case()
        before = copy.deepcopy(original)
        cfg, meta = RUNNER.resolve(original)
        self.assertEqual(original, before)
        self.assertEqual(meta['interaction'], cfg['interaction'])
        values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
        self.assertEqual(values['surface_adhesion'], 1)
        self.assertEqual(values['adhesion_work'], .0219)
        self.assertEqual(values['adhesion_range'], 6.7e-10)
        self.assertEqual(values['curvature_switch_gap'], 5e-9)
        self.assertEqual(values['curvature_cutoff_gap'], 2e-8)
        self.assertEqual(values['tangential_stiffness'], 80)
        self.assertFalse(any(key.startswith('local_') for key in values))
        self.assertEqual(RUNNER.resolve(cfg)[0], cfg)

    def test_opt_in_is_required_and_has_declared_defaults(self):
        cfg = case()
        for key in RUNNER.SURFACE_ADHESION_DEFAULTS:
            cfg['interaction'].pop(key)
        enabled, _ = RUNNER.resolve(cfg)
        for key, expected in RUNNER.SURFACE_ADHESION_DEFAULTS.items():
            self.assertEqual(enabled['interaction'][key], expected)
        cfg['interaction'].pop('surface_adhesion')
        legacy, _ = RUNNER.resolve(cfg)
        self.assertFalse(legacy['interaction']['surface_adhesion'])
        self.assertEqual(legacy['interaction']['local_gap_fraction'], 0)

    def test_invalid_material_or_mixed_inputs_rejected(self):
        for name, cfg in invalid_cases():
            with self.subTest(name=name), self.assertRaises(ValueError):
                RUNNER.resolve(cfg)

    def test_build_requires_actual_surface_capability(self):
        cfg, _ = RUNNER.resolve(case())
        for info in ({}, {'local_gap_adhesion': True}, {'surface_adhesion_version': 2}):
            with self.subTest(info=info), self.assertRaisesRegex(ValueError, 'Rebuild'):
                RUNNER.require_local_adhesion_build(cfg, info)
        RUNNER.require_local_adhesion_build(cfg, {'surface_adhesion_version': 1})


@unittest.skipUnless(shutil.which('g++'), 'C++ integration requires g++')
class CppSurfaceAdhesionInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='surface-adhesion-config-')
        cls.directory = Path(cls.temporary.name)
        source = cls.directory/'parse.cpp'
        source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  try {
    const auto c=slurry::gr_re2::parseConfig(argc,argv);
    std::cout<<std::setprecision(17)<<c.surface_adhesion<<" "<<c.adhesion_work<<" "
      <<c.adhesion_range<<" "<<c.curvature_switch_gap<<" "<<c.curvature_cutoff_gap
      <<" "<<c.roughness_gap<<" "<<c.tangential_stiffness;
  } catch(const std::exception& e) { std::cerr<<e.what();return 2; }
}
''')
        cls.executable = cls.directory/'parse'
        subprocess.run(['g++', '-std=c++17', '-I'+str(ROOT/'olb-1.9r0/src/slurry/gr_re2'),
                        str(source), '-o', str(cls.executable)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def parse_values(self, changes=None):
        cfg, _ = RUNNER.resolve(case())
        values = RUNNER.solver_values(cfg, self.directory, self.directory/'particles.csv', 0)
        values.update(changes or {})
        path = self.directory/'run.cfg'
        path.write_text(''.join('{}={}\n'.format(key,value) for key,value in values.items()))
        return subprocess.run([str(self.executable),'--config',str(path)],
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

    def test_generated_parameters_reach_cpp(self):
        result = self.parse_values()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([float(value) for value in result.stdout.split()],
                         [1,.0219,6.7e-10,5e-9,2e-8,2e-9,80])

    def test_cpp_rejects_direct_invalid_or_mixed_inputs(self):
        for changes in ({'surface_adhesion':2}, {'surface_adhesion':0}, {'adhesion_work':1e-5},
                        {'adhesion_work':'nan'}, {'adhesion_range':0}, {'adhesion_range':'inf'},
                        {'adhesion_range':4e-9}, {'curvature_switch_gap':2e-9},
                        {'curvature_cutoff_gap':5e-7}, {'rough_contact_enabled':0},
                        {'sigma_lj':4e-9}, {'local_gap_fraction':0}, {'local_gap':3e-10},
                        {'local_switch_excess_gap':2e-9}, {'local_cutoff_excess_gap':1e-8}):
            with self.subTest(changes=changes):
                result = self.parse_values(changes)
                self.assertEqual(result.returncode, 2, result.stdout)


if __name__ == '__main__':
    unittest.main()
