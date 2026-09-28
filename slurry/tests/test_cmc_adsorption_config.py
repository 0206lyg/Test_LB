"""CMC passivation, independent concentration inputs, and restart physics."""
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
SPEC = importlib.util.spec_from_file_location('cmc_runner', ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def case():
    return json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())


def resolved(inputs=None):
    original = case()
    if inputs is not None:
        original['cmc'] = inputs
    cfg, meta = RUNNER.resolve(original)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    return cfg, meta, values


def saved_checkpoint(cfg, meta):
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    values['interaction_model_version'] = RUNNER.SURFACE_ADHESION_VERSION
    immutable = {key: values[key] for key in RUNNER.SURFACE_CHECKPOINT_KEYS}
    immutable.update({key: values[key] for key in RUNNER.CMC_CHECKPOINT_KEYS if key in values})
    immutable.update(dt_s=meta['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    return {'step': 1, 'immutable_config': immutable}


class CmcInputTests(unittest.TestCase):
    def test_absent_and_zero_cmc_preserve_exact_solver_physics(self):
        pure, pure_meta, pure_values = resolved()
        self.assertNotIn('cmc', pure)
        self.assertNotIn('cmc', pure_meta)
        for inputs in ({}, {'adsorbed_g_L': 0, 'free_g_L': 0, 'q_sat': .33}):
            with self.subTest(inputs=inputs):
                cfg, meta, values = resolved(inputs)
                self.assertEqual(values, pure_values)
                self.assertEqual(cfg['interaction'], pure['interaction'])
                self.assertEqual(meta['cmc']['q'], 1.0)
                self.assertEqual(meta['interaction']['adhesion_work_J_m2'], .0219)

    def test_partial_adsorption_changes_only_extra_adhesion(self):
        _, _, pure_values = resolved()
        cfg, meta, values = resolved({'adsorbed_g_L': 1.0, 'adsorbed_saturation_g_L': 2.0})
        state = meta['cmc']
        expected_q = (1-(1-math.sqrt(.33))*.5)**2
        h0 = cfg['rough_contact']['roughness_gap_m']
        interaction = cfg['interaction']
        expected_bg = interaction['hamaker_J']/(12*math.pi*h0*h0) * (
            1-(interaction['sigma_lj_m']/h0)**6/30)
        expected_work = expected_bg+expected_q*(.0219-expected_bg)
        self.assertEqual(state['theta'], .5)
        self.assertAlmostEqual(state['q'], expected_q)
        self.assertEqual(state['background_work_J_m2'], expected_bg)
        self.assertEqual(values['adhesion_work'], expected_work)
        self.assertEqual(meta['interaction']['adhesion_work_J_m2'], expected_work)
        self.assertEqual(meta['interaction']['bare_adhesion_work_J_m2'], .0219)
        self.assertEqual(cfg['interaction']['adhesion_work_J_m2'], .0219)
        self.assertEqual({key for key in values if values[key] != pure_values.get(key)},
                         {'adhesion_work', 'cmc_contact_gap', 'cmc_cohesion_retention', 'cmc_contact_version'})
        self.assertEqual(values['cmc_contact_gap'], h0 + 1e-9)
        self.assertEqual(values['roughness_gap'], h0)

    def test_saturation_clamps_without_moving_excess_to_free_cmc(self):
        cfg, meta, values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 0})
        state = meta['cmc']
        reference = .0037*997*.44/.56
        self.assertEqual(state['adsorbed_saturation_g_L'], reference)
        self.assertEqual(state['effective_adsorbed_g_L'], reference)
        self.assertEqual(state['adsorbed_g_L'], 3.8)
        self.assertEqual(cfg['cmc']['adsorbed_g_L'], 3.8)
        self.assertTrue(state['adsorbed_clamped'])
        self.assertEqual(state['free_g_L'], 0)
        self.assertEqual(state['theta'], 1.0)
        self.assertEqual(state['q'], .33)
        self.assertEqual(values['adhesion_work'], state['effective_adhesion_work_J_m2'])

    def test_disabled_free_cohesion_records_cmc_without_physical_effect(self):
        cfg0, _, values0 = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 0})
        cfg1, meta1, values1 = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 123.0,
                                      'free_cohesion': {'enabled': False}})
        self.assertEqual(values0, values1)
        self.assertEqual(cfg0['fluid'], cfg1['fluid'])
        self.assertEqual(meta1['cmc']['free_g_L'], 123.0)
        self.assertFalse(meta1['cmc']['free_cmc_physics_enabled'])
        self.assertEqual(resolved({'free_g_L': 123.0,
                                  'free_cohesion': {'enabled': False}})[2], resolved()[2])

    def test_resolve_roundtrip_is_idempotent_and_does_not_mutate_input(self):
        original = case()
        original['cmc'] = {'adsorbed_g_L': 1.5, 'free_g_L': 2.0, '_units': 'g/L liquid'}
        before = copy.deepcopy(original)
        cfg1, meta1 = RUNNER.resolve(original)
        cfg2, meta2 = RUNNER.resolve(json.loads(json.dumps(cfg1)))
        self.assertEqual(original, before)
        self.assertEqual(cfg1, cfg2)
        self.assertEqual(meta1, meta2)
        self.assertEqual(cfg2['cmc']['_units'], 'g/L liquid')
        self.assertNotIn('_units', meta2['cmc'])
        self.assertEqual(cfg2['interaction']['adhesion_work_J_m2'], .0219)

    def test_q_endpoints_keep_background_and_exact_bare_work(self):
        cfg, meta, values = resolved({'adsorbed_g_L': 3.8, 'q_sat': 0.0})
        self.assertEqual(meta['cmc']['q'], 0.0)
        self.assertEqual(values['adhesion_work'], meta['cmc']['background_work_J_m2'])
        self.assertGreater(values['adhesion_work'], 0)
        for inputs in ({'adsorbed_g_L': 3.8, 'q_sat': 1.0, 'contact_offset_at_saturation_m': 0},
                       {'adsorbed_g_L': 0.0, 'q_sat': 0.0}):
            with self.subTest(inputs=inputs):
                self.assertEqual(resolved(inputs)[2], resolved()[2])

    def test_reference_saturation_does_not_follow_target_loading_or_fluid_density(self):
        original = case()
        original['cmc'] = {'adsorbed_g_L': 1.0}
        _, meta0 = RUNNER.resolve(original)
        original['particles']['target_solid_mass_fraction'] = .2
        original['fluid']['density_kg_m3'] = 1100.0
        _, meta1 = RUNNER.resolve(original)
        self.assertEqual(meta0['cmc'], meta1['cmc'])

    def test_invalid_cmc_inputs_are_rejected(self):
        for value in (None, [], True, 3.8, 'CMC'):
            original = case()
            original['cmc'] = value
            with self.subTest(section=value), self.assertRaises(ValueError):
                RUNNER.resolve(original)
        for key in ('adsorbed_g_L', 'free_g_L', 'q_sat', 'adsorbed_saturation_g_L',
                    'contact_offset_at_saturation_m'):
            bad_values = [-1, math.nan, math.inf, -math.inf, True, '0.33', None, []]
            if key == 'q_sat':
                bad_values.append(1.01)
            if key == 'adsorbed_saturation_g_L':
                bad_values.append(0)
            for value in bad_values:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    resolved({key: value})
        with self.assertRaisesRegex(ValueError, 'Unknown cmc'):
            resolved({'adsrobed_g_L': 1.0})

    def test_positive_adsorption_requires_surface_adhesion(self):
        legacy = case()
        legacy['interaction'].pop('surface_adhesion')
        for key in RUNNER.SURFACE_ADHESION_DEFAULTS:
            legacy['interaction'].pop(key)
        baseline, _ = RUNNER.resolve(legacy)
        legacy['cmc'] = {'adsorbed_g_L': 0, 'free_g_L': 2.0,
                         'free_cohesion': {'enabled': False}}
        cfg, meta = RUNNER.resolve(legacy)
        self.assertEqual(RUNNER.solver_values(cfg, Path('r'), Path('p'), 0),
                         RUNNER.solver_values(baseline, Path('r'), Path('p'), 0))
        self.assertIsNone(meta['cmc']['effective_adhesion_work_J_m2'])
        legacy['cmc']['adsorbed_g_L'] = 1.0
        with self.assertRaisesRegex(ValueError, 'surface_adhesion=true'):
            RUNNER.resolve(legacy)

    def test_restart_uses_effective_work_and_ignores_record_only_free_cmc(self):
        cfg, meta, _ = resolved({'adsorbed_g_L': 1.0})
        checkpoint = saved_checkpoint(cfg, meta)
        for inputs in ({'adsorbed_g_L': 1.1}, {'adsorbed_g_L': 1.0, 'q_sat': .4}):
            new_cfg, new_meta, _ = resolved(inputs)
            with self.subTest(inputs=inputs), self.assertRaisesRegex(ValueError, 'adhesion_work'):
                RUNNER.validate_restart(new_cfg, new_meta, checkpoint, 1)
        new_cfg, new_meta, _ = resolved({'adsorbed_g_L': 1.0, 'free_g_L': 16,
                                       'free_cohesion': {'enabled': False}})
        self.assertTrue(RUNNER.validate_restart(new_cfg, new_meta, checkpoint, 1))
        saturated, sat_meta, _ = resolved({'adsorbed_g_L': 3.8})
        new_cfg, new_meta, _ = resolved({'adsorbed_g_L': 9.0})
        self.assertTrue(RUNNER.validate_restart(new_cfg, new_meta, saved_checkpoint(saturated, sat_meta), 1))

    def test_coated_contact_gap_is_independent_of_free_cmc(self):
        cfg, meta, values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 0})
        self.assertEqual(values['cmc_contact_gap'], 4e-9)
        self.assertEqual(cfg['particles']['minimum_gap_m'], 4e-9)
        self.assertEqual(meta['requested_minimum_gap_m'], 3e-9)
        self.assertEqual(meta['effective_minimum_gap_m'], 4e-9)
        _, _, free_values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 13.2})
        for key in ('cmc_contact_gap', 'roughness_gap', 'sliding_friction',
                    'tangential_stiffness', 'rolling_length', 'adhesion_work'):
            self.assertEqual(values[key], free_values[key])
        cfg2, meta2 = RUNNER.resolve(json.loads(json.dumps(cfg)))
        self.assertEqual(cfg, cfg2)
        self.assertEqual(meta, meta2)

    def test_contact_plane_cannot_push_adhesion_into_background_blend(self):
        with self.assertRaisesRegex(ValueError, 'contact_gap_m.*curvature_switch'):
            resolved({'adsorbed_g_L': 3.8, 'contact_offset_at_saturation_m': 3e-9})

    def test_shipped_cmc_case_accepts_documentation_fields(self):
        original = json.loads((ROOT/'slurry/cases/gr_CMC.json').read_text())
        cfg, meta = RUNNER.resolve(original)
        self.assertEqual(cfg['cmc']['adsorbed_g_L'], 3.8)
        self.assertEqual(cfg['cmc']['free_g_L'], 0)
        self.assertEqual(meta['cmc']['q'], .33)
        self.assertTrue(meta['cmc']['adsorbed_clamped'])


class CmcRestartSelectionTests(unittest.TestCase):
    def fixture(self, run):
        checkpoint = run/'checkpoints/checkpoint_00000000000000000008'
        checkpoint.mkdir(parents=True)
        files = {'state.bin': 1, 'initial_particles.csv': 1, 'lattice_rank_0.bin': 1}
        for name in files:
            (checkpoint/name).write_bytes(b'x')
        RUNNER.write_json(checkpoint/'checkpoint.json', {
            'engine': 'pure_gr', 'format_version': 1, 'complete': True, 'step': 8,
            'ranks': 1, 'dt_s': 1e-7, 'shear_rate_s_inv': 100, 'files': files})
        return checkpoint

    def test_engine_selection_does_not_mix_case_subdirectories(self):
        with tempfile.TemporaryDirectory(prefix='cmc-restarts-') as directory:
            root = Path(directory)
            pure = self.fixture(root/'pure_gr/g000_100')
            with self.assertRaisesRegex(ValueError, 'No gr_cmc checkpoint'):
                RUNNER.restart_targets(root, 'gr_cmc')
            cmc = self.fixture(root/'gr_cmc/g000_100')
            self.assertEqual([p for p, _ in RUNNER.restart_targets(root)], [pure])
            self.assertEqual([p for p, _ in RUNNER.restart_targets(root, 'gr_cmc')], [cmc])
            self.assertEqual(RUNNER.restart_targets(cmc, 'gr_cmc')[0][0], cmc)
            with self.assertRaisesRegex(ValueError, 'Requested gr_cmc'):
                RUNNER.restart_targets(root/'pure_gr', 'gr_cmc')


@unittest.skipUnless(shutil.which('g++'), 'C++ parser integration requires g++')
class CmcCppWiringTests(unittest.TestCase):
    def test_effective_work_reaches_existing_cpp_parser(self):
        with tempfile.TemporaryDirectory(prefix='cmc-parser-') as directory:
            directory = Path(directory)
            source = directory/'parse.cpp'
            source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  const auto c=slurry::gr_re2::parseConfig(argc,argv);
  std::cout<<std::setprecision(17)<<c.adhesion_work<<" "<<c.sliding_friction;
}
''')
            executable = directory/'parse'
            subprocess.run(['g++', '-std=c++17', '-I'+str(ROOT/'olb-1.9r0/src/slurry/gr_re2'),
                            str(source), '-o', str(executable)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for inputs in ({'adsorbed_g_L': 3.8}, {'adsorbed_g_L': 3.8, 'q_sat': 0.0}):
                with self.subTest(inputs=inputs):
                    cfg, meta, values = resolved(inputs)
                    path = directory/'run.cfg'
                    path.write_text(''.join('{}={}\n'.format(k,v) for k,v in values.items()))
                    result = subprocess.run([str(executable), '--config', str(path)], check=True,
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            universal_newlines=True)
                    work, friction = map(float, result.stdout.split())
                    self.assertEqual(work, meta['cmc']['effective_adhesion_work_J_m2'])
                    self.assertEqual(friction, cfg['rough_contact']['sliding_friction'])


if __name__ == '__main__':
    unittest.main()
