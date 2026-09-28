"""Coated-contact cohesion screening, migration, and solver-policy wiring."""
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
SPEC = importlib.util.spec_from_file_location(
    'free_cmc_runner', ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
GAP = 'cmc_contact_gap'
RETENTION = 'cmc_cohesion_retention'
VERSION = 'cmc_contact_version'


def case(cmc=None):
    value = json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())
    if cmc is not None:
        value['cmc'] = cmc
    return value


def resolved(cmc=None, source=None):
    cfg, meta = RUNNER.resolve(case(cmc) if source is None else source)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    return cfg, meta, values


def checkpoint(cfg, meta, values):
    immutable = {key: values[key] for key in RUNNER.SURFACE_CHECKPOINT_KEYS
                 + RUNNER.CMC_CHECKPOINT_KEYS if key in values}
    immutable['interaction_model_version'] = RUNNER.SURFACE_ADHESION_VERSION
    immutable.update(dt_s=meta['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    return {'step': 1, 'immutable_config': immutable}


class FreeCmcConcentrationTests(unittest.TestCase):
    def test_independent_uncapped_concentration_and_exact_endpoints(self):
        for free, expected in ((0, 1), (6.6, .5), (13.2, 0), (123, 0), (1e308, 0)):
            cfg, meta, values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': free})
            state = meta['cmc']['free_cohesion']
            self.assertEqual(state['free_concentration_kg_m3'], free)
            self.assertEqual(cfg['cmc']['free_g_L'], free)
            self.assertEqual(values[RETENTION], expected)
            self.assertEqual(state['cohesion_retention'], expected)
            self.assertEqual(values[GAP], 4e-9)
            self.assertFalse(any(key.startswith('free_cmc_repulsion') for key in values))

    def test_only_screening_changes_with_free_concentration(self):
        _, _, baseline = resolved({'adsorbed_g_L': 1.2})
        _, _, values = resolved({'adsorbed_g_L': 1.2, 'free_g_L': 10})
        self.assertEqual({key for key in values if values[key] != baseline[key]}, {RETENTION})

    def test_disabled_screening_preserves_adsorbed_contact_model(self):
        _, _, baseline = resolved({'adsorbed_g_L': 1.2})
        _, meta, values = resolved({'adsorbed_g_L': 1.2, 'free_g_L': 123,
                                    'free_cohesion': {'enabled': False}})
        self.assertEqual(values, baseline)
        self.assertFalse(meta['cmc']['free_cmc_physics_enabled'])
        self.assertEqual(meta['cmc']['free_g_L'], 123)

    def test_smooth_fraction_and_endpoint_scale(self):
        source = case({'free_g_L': 3, 'free_cohesion': {'nonadhesive_concentration_g_L': 12}})
        _, _, values = resolved(source=source)
        u = .25
        self.assertAlmostEqual(values[RETENTION], 1-10*u**3+15*u**4-6*u**5)
        source['fluid']['temperature_K'] *= 1.5
        source['fluid']['density_kg_m3'] *= 1.2
        self.assertEqual(resolved(source=source)[2][RETENTION], values[RETENTION])

    def test_roundtrip_comments_and_requested_minimum_gap(self):
        source = case({'adsorbed_g_L': 3.8, 'free_g_L': 2,
                       'free_cohesion': {'_note': 'effective closure'}})
        before = copy.deepcopy(source)
        cfg, meta, _ = resolved(source=source)
        cfg2, meta2, _ = resolved(source=json.loads(json.dumps(cfg)))
        self.assertEqual(source, before)
        self.assertEqual(cfg, cfg2)
        self.assertEqual(meta, meta2)
        self.assertEqual(meta2['requested_minimum_gap_m'], 3e-9)
        self.assertEqual(meta2['effective_minimum_gap_m'], 4e-9)
        self.assertEqual(cfg['cmc']['free_cohesion']['_note'], 'effective closure')
        self.assertNotIn('_note', meta['cmc']['free_cohesion'])
        cfg['particles']['minimum_gap_m'] = 5e-9
        cfg3, meta3, _ = resolved(source=cfg)
        self.assertEqual(meta3['requested_minimum_gap_m'], 5e-9)
        self.assertEqual(cfg3['particles']['minimum_gap_m'], 5e-9)

    def test_invalid_nested_inputs_and_legacy_law_are_rejected(self):
        for value in (None, [], True, 1, 'yes'):
            with self.subTest(section=value), self.assertRaises(ValueError):
                resolved({'free_cohesion': value})
        for value in (0, 1, 'true', None):
            with self.subTest(enabled=value), self.assertRaises(ValueError):
                resolved({'free_cohesion': {'enabled': value}})
        for value in (0, -1, math.nan, math.inf, True, '13.2', None, []):
            with self.subTest(endpoint=value), self.assertRaises(ValueError):
                resolved({'free_cohesion': {'nonadhesive_concentration_g_L': value}})
        with self.assertRaisesRegex(ValueError, 'Unknown cmc.free_cohesion'):
            resolved({'free_cohesion': {'strength': 1}})
        for old in ({}, {'enabled': False}, {'strength': 1}):
            with self.subTest(old=old), self.assertRaisesRegex(ValueError, 'obsolete exponential'):
                resolved({'free_repulsion': old})


class FreeCmcCompatibilityTests(unittest.TestCase):
    def test_build_capabilities(self):
        cfg, _, _ = resolved({'adsorbed_g_L': 1})
        for info in ({}, {'surface_adhesion_version': 1, 'pass_max_version': 1},
                     {'surface_adhesion_version': 1, 'pass_max_version': 1, VERSION: 2}):
            with self.subTest(info=info), self.assertRaisesRegex(ValueError, 'Rebuild'):
                RUNNER.require_local_adhesion_build(cfg, info)
        RUNNER.require_local_adhesion_build(cfg, {
            'surface_adhesion_version': 1, 'pass_max_version': 1, VERSION: 1})
        zero, _, _ = resolved({'free_g_L': 0})
        RUNNER.require_local_adhesion_build(zero, {'surface_adhesion_version': 1, 'pass_max_version': 1})

    def test_same_physics_and_old_pure_checkpoints_restart(self):
        for cmc in (None, {'free_g_L': 0}, {'adsorbed_g_L': 3.8, 'free_g_L': 13.2}):
            with self.subTest(cmc=cmc):
                cfg, meta, values = resolved(cmc)
                saved = checkpoint(cfg, meta, values)
                self.assertTrue(RUNNER.validate_restart(cfg, meta, saved, 1))
                cfg['numerics']['pass_max'] = 0
                self.assertTrue(RUNNER.validate_restart(cfg, meta, saved, 1))

    def test_changed_cohesion_contact_or_old_exponential_cannot_restart(self):
        cfg, meta, values = resolved({'adsorbed_g_L': 1, 'free_g_L': 2})
        saved = checkpoint(cfg, meta, values)
        for cmc in ({'adsorbed_g_L': 1, 'free_g_L': 0},
                    {'adsorbed_g_L': 1, 'free_g_L': 3},
                    {'adsorbed_g_L': 1, 'free_g_L': 2, 'contact_offset_at_saturation_m': 1e-9},
                    {'adsorbed_g_L': 1, 'free_g_L': 2, 'free_cohesion': {'enabled': False}}):
            new_cfg, new_meta, _ = resolved(cmc)
            with self.subTest(cmc=cmc), self.assertRaises(ValueError):
                RUNNER.validate_restart(new_cfg, new_meta, saved, 1)
        saved['immutable_config']['free_cmc_repulsion_pressure'] = 0
        with self.assertRaisesRegex(ValueError, 'obsolete'):
            RUNNER.validate_restart(cfg, meta, saved, 1)

    def test_incomplete_contact_checkpoint_rejected(self):
        cfg, meta, values = resolved({'adsorbed_g_L': 1, 'free_g_L': 2})
        saved = checkpoint(cfg, meta, values)
        for key in RUNNER.CMC_CHECKPOINT_KEYS:
            incomplete = copy.deepcopy(saved)
            del incomplete['immutable_config'][key]
            with self.subTest(key=key), self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, meta, incomplete, 1)

    def test_pass_max_defaults_and_rejects_nonbinary_values(self):
        cfg, meta, values = resolved()
        self.assertEqual(cfg['numerics']['pass_max'], 1)
        self.assertEqual(meta['pass_max'], 1)
        self.assertEqual(values['pass_max'], 1)
        source = case()
        source['numerics']['pass_max'] = 0
        self.assertEqual(resolved(source=source)[2]['pass_max'], 0)
        for value in (-1, 2, True, 1.0, '1', None):
            source['numerics']['pass_max'] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'pass_max'):
                resolved(source=source)


@unittest.skipUnless(shutil.which('g++'), 'C++ parser integration requires g++')
class FreeCmcCppWiringTests(unittest.TestCase):
    def test_generated_contact_cohesion_and_pass_max_reach_cpp(self):
        with tempfile.TemporaryDirectory(prefix='free-cmc-config-') as temporary:
            directory = Path(temporary)
            source = directory/'parse.cpp'
            source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  const auto c=slurry::gr_re2::parseConfig(argc,argv);
  std::cout<<std::setprecision(17)<<c.cmc_contact_gap<<" "
    <<c.cmc_cohesion_retention<<" "<<c.cmc_contact_version<<" "
    <<c.pass_max<<" "<<c.adhesion_work;
}
''')
            executable = directory/'parse'
            subprocess.run(['g++', '-std=c++17', '-I'+str(ROOT/'olb-1.9r0/src/slurry/gr_re2'),
                            str(source), '-o', str(executable)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for free in (0, 13.2):
                cfg, meta, values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': free})
                path = directory/'run.cfg'
                path.write_text(''.join('{}={}\n'.format(key, value) for key, value in values.items()))
                result = subprocess.run([str(executable), '--config', str(path)], check=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                gap, retention, version, pass_max, work = map(float, result.stdout.split())
                self.assertEqual((gap, retention, version, pass_max),
                                 (values[GAP], values[RETENTION], 1, 1))
                self.assertEqual(work, meta['cmc']['effective_adhesion_work_J_m2'])


if __name__ == '__main__':
    unittest.main()
