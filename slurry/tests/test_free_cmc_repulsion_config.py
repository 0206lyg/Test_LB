"""Free-CMC SI conversion, solver wiring, and restart compatibility."""
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
VERSION = 'free_cmc_repulsion_version'
R = 8.31446261815324


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
                 if key in values}
    immutable['interaction_model_version'] = RUNNER.SURFACE_ADHESION_VERSION
    immutable.update(dt_s=meta['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    if PRESSURE in values:
        immutable.update({PRESSURE: values[PRESSURE], LENGTH: values[LENGTH], VERSION: 1})
    return {'step': 1, 'immutable_config': immutable}


class FreeCmcConcentrationTests(unittest.TestCase):
    def test_g_per_litre_is_kg_per_cubic_metre_without_extra_thousand(self):
        cfg, meta, values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 13.2})
        state = meta['cmc']['free_repulsion']
        expected = .5*.7/.218*R*298.15*13.2
        self.assertEqual(state['free_concentration_kg_m3'], 13.2)
        self.assertAlmostEqual(state['repeat_unit_concentration_mol_m3'], 13.2/.218)
        self.assertAlmostEqual(state['nominal_charge_concentration_mol_m3'], .7*13.2/.218)
        self.assertAlmostEqual(state['bulk_osmotic_pressure_Pa'], expected)
        self.assertAlmostEqual(state['effective_repulsion_pressure_Pa'], expected)
        self.assertAlmostEqual(values[PRESSURE], expected)
        self.assertEqual(values[LENGTH], 5e-9)
        self.assertAlmostEqual(state['surface_energy_at_contact_J_m2']/(expected*5e-9), 1)
        self.assertAlmostEqual(state['pair_energy_per_derjaguin_length_J_m'] /
                               (math.pi*expected*(5e-9)**2), 1)
        self.assertEqual(state['model_version'], 1)
        self.assertEqual(state['temperature_K'], cfg['fluid']['temperature_K'])
        self.assertTrue(state['active'])
        self.assertTrue(meta['cmc']['free_cmc_physics_enabled'])

    def test_free_concentration_is_independent_and_not_capped_by_adsorption(self):
        original = case({'adsorbed_g_L': 0, 'free_g_L': 123.0})
        _, baseline, _ = resolved(source=original)
        original['cmc'].update(adsorbed_g_L=99, adsorbed_saturation_g_L=2)
        original['particles']['target_solid_mass_fraction'] = .2
        original['fluid']['density_kg_m3'] = 1100
        _, changed, _ = resolved(source=original)
        self.assertEqual(changed['cmc']['free_g_L'], 123)
        self.assertEqual(changed['cmc']['free_repulsion'], baseline['cmc']['free_repulsion'])

    def test_pressure_scales_with_concentration_strength_and_temperature(self):
        _, _, reference = resolved({'free_g_L': 1.2})
        source = case({'free_g_L': 2.4, 'free_repulsion': {'strength': 3}})
        source['fluid']['temperature_K'] *= 1.5
        _, _, values = resolved(source=source)
        self.assertAlmostEqual(values[PRESSURE]/reference[PRESSURE], 9)

    def test_length_changes_energy_without_changing_osmotic_pressure(self):
        _, m1, v1 = resolved({'free_g_L': 2})
        _, m2, v2 = resolved({'free_g_L': 2, 'free_repulsion': {'decay_length_m': 1e-8}})
        self.assertEqual(v1[PRESSURE], v2[PRESSURE])
        a, b = m1['cmc']['free_repulsion'], m2['cmc']['free_repulsion']
        self.assertAlmostEqual(b['surface_energy_at_contact_J_m2']/a['surface_energy_at_contact_J_m2'], 2)
        self.assertAlmostEqual(b['pair_energy_per_derjaguin_length_J_m']/a['pair_energy_per_derjaguin_length_J_m'], 4)

    def test_zero_disabled_and_zero_strength_preserve_solver_dictionary(self):
        _, _, baseline = resolved({'adsorbed_g_L': 1.2})
        for cmc in ({'free_g_L': 0},
                    {'free_g_L': 10, 'free_repulsion': {'enabled': False}},
                    {'free_g_L': 10, 'free_repulsion': {'strength': 0}}):
            with self.subTest(cmc=cmc):
                _, meta, values = resolved(dict(cmc, adsorbed_g_L=1.2))
                self.assertEqual(values, baseline)
                self.assertFalse(meta['cmc']['free_cmc_physics_enabled'])
                self.assertFalse(meta['cmc']['free_repulsion']['active'])

    def test_active_repulsion_does_not_change_adhesion_or_contact_parameters(self):
        _, _, baseline = resolved({'adsorbed_g_L': 1.2})
        _, _, values = resolved({'adsorbed_g_L': 1.2, 'free_g_L': 10})
        self.assertEqual({key: value for key, value in values.items()
                          if key not in (PRESSURE, LENGTH)}, baseline)

    def test_roundtrip_and_comment_fields_preserve_user_input(self):
        source = case({'free_g_L': 2, 'free_repulsion': {'_note': 'trial barrier', 'strength': .7}})
        before = copy.deepcopy(source)
        cfg, meta, _ = resolved(source=source)
        cfg2, meta2, _ = resolved(source=json.loads(json.dumps(cfg)))
        self.assertEqual(source, before)
        self.assertEqual(cfg, cfg2)
        self.assertEqual(meta, meta2)
        self.assertEqual(cfg['cmc']['free_repulsion']['_note'], 'trial barrier')
        self.assertNotIn('_note', meta['cmc']['free_repulsion'])

    def test_invalid_nested_inputs_are_rejected(self):
        for value in (None, [], True, 1, 'yes'):
            with self.subTest(section=value), self.assertRaises(ValueError):
                resolved({'free_g_L': 1, 'free_repulsion': value})
        for value in (0, 1, 'true', None):
            with self.subTest(enabled=value), self.assertRaises(ValueError):
                resolved({'free_repulsion': {'enabled': value}})
        for key in ('strength', 'decay_length_m', 'degree_of_substitution',
                    'repeat_unit_molar_mass_kg_mol', 'osmotic_coefficient'):
            for value in (-1, math.nan, math.inf, -math.inf, True, '0.5', None, []):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    resolved({'free_g_L': 1, 'free_repulsion': {key: value}})
        for key in ('decay_length_m', 'repeat_unit_molar_mass_kg_mol'):
            with self.subTest(zero=key), self.assertRaises(ValueError):
                resolved({'free_g_L': 1, 'free_repulsion': {key: 0}})
        with self.assertRaisesRegex(ValueError, 'Unknown cmc.free_repulsion'):
            resolved({'free_repulsion': {'strenght': 1}})

    def test_nonfinite_derived_pressure_is_rejected(self):
        for cmc in ({'free_g_L': 1e308},
                    {'free_g_L': 1, 'free_repulsion': {'strength': 1e308}},
                    {'free_g_L': 1, 'free_repulsion': {'repeat_unit_molar_mass_kg_mol': 1e-320}}):
            with self.subTest(cmc=cmc), self.assertRaises(ValueError):
                resolved(cmc)

    def test_active_energy_underflow_is_rejected_instead_of_silently_disabling(self):
        with self.assertRaises(ValueError):
            resolved({'free_g_L': 1e-320})

    def test_repulsion_requires_the_contact_gap_reference(self):
        source = case({'free_g_L': 2})
        source['rough_contact']['enabled'] = False
        source['interaction'].pop('surface_adhesion')
        for key in RUNNER.SURFACE_ADHESION_DEFAULTS:
            source['interaction'].pop(key)
        with self.assertRaises(ValueError):
            resolved(source=source)


class FreeCmcCompatibilityTests(unittest.TestCase):
    def test_active_case_requires_new_executable_capability(self):
        cfg, _, _ = resolved({'free_g_L': 1})
        for version in (None, 0, 2):
            info = {'surface_adhesion_version': 1}
            if version is not None:
                info[VERSION] = version
            with self.subTest(version=version), self.assertRaisesRegex(ValueError, 'Rebuild'):
                RUNNER.require_local_adhesion_build(cfg, info)
        RUNNER.require_local_adhesion_build(cfg, {'surface_adhesion_version': 1, VERSION: 1})
        zero, _, _ = resolved({'free_g_L': 0})
        RUNNER.require_local_adhesion_build(zero, {'surface_adhesion_version': 1})

    def test_unchanged_active_and_old_zero_checkpoints_can_restart(self):
        for cmc in (None, {'free_g_L': 0}, {'adsorbed_g_L': 3.8, 'free_g_L': 13.2}):
            with self.subTest(cmc=cmc):
                cfg, meta, values = resolved(cmc)
                self.assertTrue(RUNNER.validate_restart(cfg, meta, checkpoint(cfg, meta, values), 1))

    def test_changed_repulsion_is_rejected_including_enable_and_disable(self):
        cfg, meta, values = resolved({'free_g_L': 2})
        saved = checkpoint(cfg, meta, values)
        for cmc in ({'free_g_L': 0}, {'free_g_L': 3},
                    {'free_g_L': 2, 'free_repulsion': {'enabled': False}},
                    {'free_g_L': 2, 'free_repulsion': {'strength': .5}},
                    {'free_g_L': 2, 'free_repulsion': {'decay_length_m': 1e-8}}):
            new_cfg, new_meta, _ = resolved(cmc)
            with self.subTest(cmc=cmc), self.assertRaises(ValueError):
                RUNNER.validate_restart(new_cfg, new_meta, saved, 1)
        zero, zero_meta, zero_values = resolved({'free_g_L': 0})
        with self.assertRaises(ValueError):
            RUNNER.validate_restart(cfg, meta, checkpoint(zero, zero_meta, zero_values), 1)

    def test_incomplete_or_wrong_version_active_checkpoint_is_rejected(self):
        cfg, meta, values = resolved({'free_g_L': 2})
        saved = checkpoint(cfg, meta, values)
        for key in (PRESSURE, LENGTH, VERSION):
            incomplete = copy.deepcopy(saved)
            del incomplete['immutable_config'][key]
            with self.subTest(missing=key), self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, meta, incomplete, 1)
        saved['immutable_config'][VERSION] = 2
        with self.assertRaises(ValueError):
            RUNNER.validate_restart(cfg, meta, saved, 1)


@unittest.skipUnless(shutil.which('g++'), 'C++ parser integration requires g++')
class FreeCmcCppWiringTests(unittest.TestCase):
    def test_generated_pressure_and_length_reach_cpp_parser(self):
        with tempfile.TemporaryDirectory(prefix='free-cmc-config-') as temporary:
            directory = Path(temporary)
            source = directory/'parse.cpp'
            source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  const auto c=slurry::gr_re2::parseConfig(argc,argv);
  std::cout<<std::setprecision(17)<<c.free_cmc_repulsion_pressure<<" "
    <<c.free_cmc_repulsion_length<<" "<<c.adhesion_work;
}
''')
            executable = directory/'parse'
            subprocess.run(['g++', '-std=c++17', '-I'+str(ROOT/'olb-1.9r0/src/slurry/gr_re2'),
                            str(source), '-o', str(executable)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            _, meta, values = resolved({'adsorbed_g_L': 3.8, 'free_g_L': 13.2,
                                       'free_repulsion': {'strength': .8, 'decay_length_m': 6e-9}})
            path = directory/'run.cfg'
            path.write_text(''.join('{}={}\n'.format(key, value) for key, value in values.items()))
            result = subprocess.run([str(executable), '--config', str(path)], check=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
            pressure, length, work = map(float, result.stdout.split())
            self.assertEqual(pressure, values[PRESSURE])
            self.assertEqual(length, values[LENGTH])
            self.assertEqual(work, meta['cmc']['effective_adhesion_work_J_m2'])


if __name__ == '__main__':
    unittest.main()
