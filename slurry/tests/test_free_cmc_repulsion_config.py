"""CMC concentration-to-force wiring, migration and checkpoint compatibility."""
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
PRESSURE = 'free_cmc_repulsion_pressure'
LENGTH = 'free_cmc_repulsion_length'
WORK = 'free_cmc_inner_repulsion_work'
RANGE = 'free_cmc_inner_repulsion_range'
POWER = 'free_cmc_inner_repulsion_power'
INNER = {'inner_work_per_g_L_J_m2': .000543218,
         'inner_range_m': 6.7e-10, 'inner_exponent': 2.1}


def case(cmc=None):
    value = json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())
    if cmc is not None:
        value['cmc'] = cmc
    return value


def resolved(cmc=None, source=None):
    cfg, meta = RUNNER.resolve(case(cmc) if source is None else source)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    return cfg, meta, values


def inputs(free=13.101578571428572, **repulsion):
    params = dict(INNER)
    params.update(repulsion)
    return {'adsorbed_g_L': 2.8984214285714285, 'free_g_L': free,
            'free_repulsion': params}


def checkpoint(cfg, meta, values):
    values = dict(values, interaction_model_version=RUNNER.SURFACE_ADHESION_VERSION)
    if PRESSURE in values:
        values['free_cmc_repulsion_version'] = RUNNER.FREE_CMC_REPULSION_VERSION
    if WORK in values:
        values['free_cmc_inner_repulsion_version'] = RUNNER.FREE_CMC_INNER_REPULSION_VERSION
    immutable = {key: values[key] for key in RUNNER.SURFACE_CHECKPOINT_KEYS
                 + RUNNER.FREE_CMC_CHECKPOINT_KEYS + RUNNER.INNER_CMC_CHECKPOINT_KEYS if key in values}
    immutable.update(dt_s=meta['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    return {'step': 1, 'immutable_config': immutable}


class FreeCmcConcentrationTests(unittest.TestCase):
    def test_independent_concentrations_and_both_unit_conversions(self):
        for free in (1.1015785714285715, 13.101578571428572, 26.2):
            cfg, meta, values = resolved(inputs(free))
            expected_pressure = .5*.7/.218*8.31446261815324*298.15*free
            state = meta['cmc']['free_repulsion']
            self.assertEqual(state['free_concentration_kg_m3'], free)
            self.assertEqual(cfg['cmc']['free_g_L'], free)
            self.assertAlmostEqual(values[PRESSURE], expected_pressure)
            self.assertEqual(values[WORK], .000543218*free)
            self.assertEqual(values[RANGE], 6.7e-10)
            self.assertEqual(values[POWER], 2.1)
            self.assertEqual(values['roughness_gap'], 2e-9)
            self.assertFalse(any(key.startswith('cmc_') for key in values))

    def test_free_concentration_changes_only_two_force_amplitudes(self):
        _, _, low = resolved(inputs(1.1))
        _, _, high = resolved(inputs(13.1))
        self.assertEqual({key for key in high if high[key] != low[key]}, {PRESSURE, WORK})

    def test_omitted_inner_work_reproduces_7004_outer_only_law(self):
        _, state, values = resolved({'adsorbed_g_L': 2.8, 'free_g_L': 13.2})
        self.assertGreater(values[PRESSURE], 0)
        self.assertNotIn(WORK, values)
        self.assertEqual(state['cmc']['free_repulsion']['inner_work_per_g_L_J_m2'], 0)
        _, _, explicit_zero = resolved(inputs(13.2, inner_work_per_g_L_J_m2=0))
        self.assertEqual(values[PRESSURE], explicit_zero[PRESSURE])

    def test_zero_and_disabled_terms_leave_adsorption_only_physics(self):
        _, _, baseline = resolved({'adsorbed_g_L': 2.8984214285714285})
        for cmc in (inputs(0), inputs(13.1, enabled=False),
                    inputs(13.1, strength=0, inner_work_per_g_L_J_m2=0)):
            cfg, meta, values = resolved(cmc)
            self.assertEqual(values, baseline)
            self.assertFalse(meta['cmc']['free_cmc_physics_enabled'])
            self.assertEqual(cfg['cmc']['free_g_L'], cmc['free_g_L'])

    def test_outer_strength_does_not_scale_inner_and_temperature_only_scales_outer(self):
        cfg, _, baseline = resolved(inputs())
        _, _, values = resolved(inputs(strength=0))
        self.assertNotIn(PRESSURE, values)
        self.assertEqual(values[WORK], baseline[WORK])
        cfg['fluid']['temperature_K'] *= 1.5
        _, _, warmer = resolved(source=cfg)
        self.assertAlmostEqual(warmer[PRESSURE], 1.5*baseline[PRESSURE])
        self.assertEqual(warmer[WORK], baseline[WORK])

    def test_resolved_config_roundtrip_preserves_explicit_placement_and_comments(self):
        source = case(inputs(_note='trial law'))
        source['particles']['minimum_gap_m'] = 1e-8
        before = copy.deepcopy(source)
        cfg, meta, _ = resolved(source=source)
        cfg2, meta2, _ = resolved(source=json.loads(json.dumps(cfg)))
        self.assertEqual(source, before)
        self.assertEqual((cfg, meta), (cfg2, meta2))
        self.assertEqual(cfg['particles']['minimum_gap_m'], 1e-8)
        self.assertEqual(cfg['cmc']['free_repulsion']['_note'], 'trial law')
        self.assertNotIn('_note', meta['cmc']['free_repulsion'])

    def test_invalid_repulsion_values_and_conflicting_legacy_model_rejected(self):
        for value in (None, [], True, 1, 'yes'):
            with self.subTest(section=value), self.assertRaises(ValueError):
                resolved({'free_repulsion': value})
        for value in (0, 1, 'true', None):
            with self.subTest(enabled=value), self.assertRaises(ValueError):
                resolved(inputs(enabled=value))
        for key in ('strength', 'decay_length_m', 'degree_of_substitution',
                    'repeat_unit_molar_mass_kg_mol', 'osmotic_coefficient',
                    'inner_work_per_g_L_J_m2', 'inner_range_m', 'inner_exponent'):
            for value in (-1, math.nan, math.inf, True, '1', None):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    resolved(inputs(**{key: value}))
        for key, value in (('decay_length_m', 0), ('inner_range_m', 0),
                           ('repeat_unit_molar_mass_kg_mol', 0),
                           ('degree_of_substitution', 3.1), ('inner_exponent', 2)):
            with self.subTest(key=key), self.assertRaises(ValueError):
                resolved(inputs(**{key: value}))
        with self.assertRaisesRegex(ValueError, 'Unknown cmc.free_repulsion'):
            resolved(inputs(inner_power=2.1))
        for old in ({'free_cohesion': {}}, {'free_cohesion': {'enabled': False}},
                    {'contact_offset_at_saturation_m': 0}):
            with self.subTest(old=old), self.assertRaisesRegex(ValueError, 'legacy cmc fields'):
                resolved(old)
        with self.assertRaises(ValueError):
            resolved(inputs(1e308))
        with self.assertRaisesRegex(ValueError, 'inner_range_m.*adhesion_range_m'):
            resolved(inputs(inner_range_m=1e-9))
        with self.assertRaises(ValueError):
            resolved(inputs(inner_work_per_g_L_J_m2=1e300))


class FreeCmcCompatibilityTests(unittest.TestCase):
    def test_build_requires_inner_capability_only_when_enabled(self):
        cfg, _, _ = resolved(inputs())
        old = {'surface_adhesion_version': 1, 'pass_max_version': 1,
               'free_cmc_repulsion_version': 1}
        with self.assertRaisesRegex(ValueError, 'inner.*Rebuild'):
            RUNNER.require_local_adhesion_build(cfg, old)
        RUNNER.require_local_adhesion_build(cfg, dict(old, free_cmc_inner_repulsion_version=1))
        cfg, _, _ = resolved(inputs(inner_work_per_g_L_J_m2=0))
        RUNNER.require_local_adhesion_build(cfg, old)
        cfg, _, _ = resolved(inputs(0))
        RUNNER.require_local_adhesion_build(cfg, {'surface_adhesion_version': 1, 'pass_max_version': 1})

    def test_same_physics_and_old_pure_or_outer_checkpoints_restart(self):
        for cmc in (None, inputs(0), inputs(), inputs(inner_work_per_g_L_J_m2=0), inputs(strength=0)):
            with self.subTest(cmc=cmc):
                cfg, meta, values = resolved(cmc)
                saved = checkpoint(cfg, meta, values)
                self.assertTrue(RUNNER.validate_restart(cfg, meta, saved, 1))
                cfg['numerics']['pass_max'] = 0
                self.assertTrue(RUNNER.validate_restart(cfg, meta, saved, 1))

    def test_changed_inner_and_outer_laws_cannot_restart(self):
        cfg, meta, values = resolved(inputs())
        saved = checkpoint(cfg, meta, values)
        for changes in ({'inner_work_per_g_L_J_m2': 0}, {'inner_work_per_g_L_J_m2': .0006},
                        {'inner_range_m': 5e-10}, {'inner_exponent': 2.2},
                        {'decay_length_m': 6e-9}, {'strength': 0}, {'enabled': False}):
            new_cfg, new_meta, _ = resolved(inputs(**changes))
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                RUNNER.validate_restart(new_cfg, new_meta, saved, 1)
        saved['immutable_config']['cmc_contact_version'] = 1
        with self.assertRaisesRegex(ValueError, 'coated-contact/free-cohesion'):
            RUNNER.validate_restart(cfg, meta, saved, 1)

    def test_incomplete_or_wrong_versions_rejected(self):
        cfg, meta, values = resolved(inputs())
        saved = checkpoint(cfg, meta, values)
        for key in RUNNER.FREE_CMC_CHECKPOINT_KEYS + RUNNER.INNER_CMC_CHECKPOINT_KEYS:
            broken = copy.deepcopy(saved)
            del broken['immutable_config'][key]
            with self.subTest(key=key), self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, meta, broken, 1)
        for key in ('free_cmc_repulsion_version', 'free_cmc_inner_repulsion_version'):
            broken = copy.deepcopy(saved)
            broken['immutable_config'][key] = 99
            with self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, meta, broken, 1)

    def test_pass_max_policy_is_preserved(self):
        cfg, meta, values = resolved()
        self.assertEqual((cfg['numerics']['pass_max'], meta['pass_max'], values['pass_max']), (1, 1, 1))
        cfg['numerics']['pass_max'] = 0
        self.assertEqual(resolved(source=cfg)[2]['pass_max'], 0)
        for value in (-1, 2, True, 1.0, '1', None):
            cfg['numerics']['pass_max'] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'pass_max'):
                resolved(source=cfg)


@unittest.skipUnless(shutil.which('g++'), 'C++ parser integration requires g++')
class FreeCmcCppWiringTests(unittest.TestCase):
    def test_derived_forces_reach_cpp_without_changing_contact_or_solver(self):
        with tempfile.TemporaryDirectory(prefix='free-cmc-config-') as temporary:
            directory = Path(temporary)
            source = directory/'parse.cpp'
            source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  const auto c=slurry::gr_re2::parseConfig(argc,argv);
  std::cout<<std::setprecision(17)<<c.free_cmc_repulsion_pressure<<" "
    <<c.free_cmc_repulsion_length<<" "<<c.free_cmc_inner_repulsion_work<<" "
    <<c.free_cmc_inner_repulsion_range<<" "<<c.free_cmc_inner_repulsion_power<<" "
    <<c.roughness_gap<<" "<<c.cmc_contact_version<<" "<<c.cmc_cohesion_retention<<" "
    <<c.pass_max<<" "<<c.adhesion_work;
}
''')
            executable = directory/'parse'
            subprocess.run(['g++', '-std=c++17', '-I'+str(ROOT/'olb-1.9r0/src/slurry/gr_re2'),
                            str(source), '-o', str(executable)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for free in (0, 1.1015785714285715, 13.101578571428572):
                cfg, meta, values = resolved(inputs(free))
                path = directory/'run.cfg'
                path.write_text(''.join('{}={}\n'.format(key, value) for key, value in values.items()))
                result = subprocess.run([str(executable), '--config', str(path)], check=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                parsed = list(map(float, result.stdout.split()))
                self.assertEqual(parsed[:5], [values.get(PRESSURE, 0), values.get(LENGTH, 5e-9),
                                              values.get(WORK, 0), values.get(RANGE, 6.7e-10),
                                              values.get(POWER, 2.1)])
                self.assertEqual(parsed[5:9], [2e-9, 0, 1, 1])
                self.assertEqual(parsed[9], meta['cmc']['effective_adhesion_work_J_m2'])


if __name__ == '__main__':
    unittest.main()
