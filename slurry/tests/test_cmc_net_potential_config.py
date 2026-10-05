"""Net CMC activation, legacy equivalence, executable gating and restart safety."""
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
    'net_cmc_runner', ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
START = 1.1015785714285715
FULL = 13.101578571428572
NET_KEYS = RUNNER.NET_CMC_CHECKPOINT_KEYS[:-1]


def case(name='gr_CMC.json'):
    return json.loads((ROOT/'slurry/cases'/name).read_text())


def resolved(source=None):
    cfg, meta = RUNNER.resolve(case() if source is None else source)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    return cfg, meta, values


def checkpoint(cfg, meta, values):
    immutable = {}
    for physical, version, keys in (
            ('surface_adhesion', 'interaction_model_version', RUNNER.SURFACE_CHECKPOINT_KEYS),
            ('free_cmc_repulsion_pressure', 'free_cmc_repulsion_version', RUNNER.FREE_CMC_CHECKPOINT_KEYS),
            ('free_cmc_inner_repulsion_work', 'free_cmc_inner_repulsion_version', RUNNER.INNER_CMC_CHECKPOINT_KEYS),
            ('cmc_net_blend', 'cmc_net_potential_version', RUNNER.NET_CMC_CHECKPOINT_KEYS)):
        if values.get(physical, 0) > 0:
            immutable.update({key: values[key] for key in keys if key != version})
            immutable[version] = 1
    immutable.update(dt_s=meta['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    return {'step': 1, 'immutable_config': immutable}


class NetCmcConfigTests(unittest.TestCase):
    def test_activation_endpoints_and_intermediate_legacy_branch_use_same_free_cmc(self):
        for free, expected in ((0, 0), (START, 0), ((START+FULL)/2, .5),
                               (FULL, 1), (FULL+12, 1)):
            source = case()
            source['cmc']['free_g_L'] = free
            _, meta, values = resolved(source)
            net = meta['cmc']['net_potential']
            with self.subTest(free=free):
                self.assertEqual(net['blend'], expected)
                self.assertEqual(net['re2_tail_weight'], 1-expected)
                self.assertEqual(values.get('cmc_net_blend', 0), expected)
                self.assertEqual(set(NET_KEYS).intersection(values), set(NET_KEYS) if expected else set())
                self.assertEqual(meta['cmc']['free_repulsion']['free_concentration_kg_m3'], free)
                if free:
                    self.assertAlmostEqual(values['free_cmc_repulsion_pressure'],
                                           .5*.7/.218*RUNNER.MOLAR_GAS_CONSTANT*298.15*free)
                    self.assertEqual(values['free_cmc_inner_repulsion_work'], .000543218*free)

    def test_shipped_4g_and_pure_gr_remain_identical_with_inactive_net_option(self):
        four = case('gr_CMC_4g_L.json')
        self.assertNotIn('net_potential', four['cmc'])
        _, _, baseline = resolved(four)
        enabled = copy.deepcopy(four)
        enabled['cmc']['net_potential'] = {'enabled': True}
        self.assertEqual(resolved(enabled)[2], baseline)
        self.assertFalse(set(NET_KEYS).intersection(baseline))
        pure = case('pure_gr.json')
        _, _, baseline = resolved(pure)
        pure['cmc'] = {'adsorbed_g_L': 0, 'free_g_L': 0,
                       'net_potential': {'enabled': True}}
        self.assertEqual(resolved(pure)[2], baseline)

    def test_disabled_net_restores_full_legacy_config_and_shipped16_changes_only_named_settings(self):
        source = case()
        source['cmc']['net_potential']['enabled'] = False
        cfg, _, disabled = resolved(source)
        del source['cmc']['net_potential']
        self.assertEqual(resolved(source)[2], disabled)
        self.assertGreater(disabled['free_cmc_repulsion_pressure'], 0)
        self.assertGreater(disabled['free_cmc_inner_repulsion_work'], 0)
        four = case('gr_CMC_4g_L.json')
        changed = {key for key, value in four['numerics'].items()
                   if not key.startswith('_') and cfg['numerics'][key] != value}
        self.assertEqual(changed, {'particle_force_absolute_tolerance_N',
                                   'particle_torque_absolute_tolerance_N_m'})
        self.assertEqual(cfg['numerics']['particle_force_absolute_tolerance_N'], 1e-15)
        self.assertEqual(cfg['numerics']['particle_torque_absolute_tolerance_N_m'], 1.65e-21)
        self.assertEqual(cfg['rough_contact'], four['rough_contact'])

    def test_reference_geometry_metadata_and_config_roundtrip(self):
        source = case()
        source['cmc']['net_potential']['_note'] = 'trial hypothesis'
        before = copy.deepcopy(source)
        cfg, meta, values = resolved(source)
        self.assertEqual(source, before)
        self.assertEqual(resolved(json.loads(json.dumps(cfg)))[:2], (cfg, meta))
        net = meta['cmc']['net_potential']
        self.assertAlmostEqual(net['reference_length_m']/13.6125e-6, 1)
        self.assertEqual(values['cmc_net_reference_length'], net['reference_length_m'])
        self.assertAlmostEqual(net['saddle_gap_m']/3e-9, 1)
        self.assertAlmostEqual(net['cutoff_gap_m']/9e-9, 1)
        self.assertEqual(net['legacy_potential_weight'], 0)
        self.assertIn('RE2 tail', net['model'])
        self.assertIn('same_free_g_L', net['interpolation'])
        self.assertNotIn('_note', net)
        cfg['particles']['thickness_m'] *= 2
        self.assertAlmostEqual(resolved(cfg)[2]['cmc_net_reference_length']/values['cmc_net_reference_length'], .5)

    def test_invalid_inputs_and_unrepresentable_scales_are_rejected(self):
        for bad in (None, [], True, 1, 'enabled'):
            source = case()
            source['cmc']['net_potential'] = bad
            with self.subTest(section=bad), self.assertRaises(ValueError):
                resolved(source)
        for key in RUNNER.CMC_NET_POTENTIAL_DEFAULTS:
            bad_values = (0, 1, 'true', None) if key == 'enabled' else (-1, math.nan, math.inf, True, None)
            for bad in bad_values:
                source = case()
                source['cmc']['net_potential'][key] = bad
                with self.subTest(key=key, value=bad), self.assertRaises(ValueError):
                    resolved(source)
        for updates in ({'unknown': 1}, {'full_free_g_L': START},
                        {'start_free_g_L': FULL+1}, {'contact_force_N': 0},
                        {'barrier_force_N': 0}, {'repulsion_range_m': 1e-6},
                        {'attraction_range_m': 1e-300}, {'contact_force_N': 1e-320},
                        {'repulsion_range_m': 1e308}):
            source = case()
            source['cmc']['net_potential'].update(updates)
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                resolved(source)
        cfg, _, _ = resolved()
        cfg['particles'].update(diameter_m=1e300, thickness_m=1e-300)
        with self.assertRaisesRegex(ValueError, 'reference_length_m'):
            RUNNER.resolve_net_cmc(cfg['cmc'], cfg)
        for section, key in (('rough_contact', 'enabled'), ('interaction', 'surface_adhesion')):
            source = case()
            source[section][key] = False
            with self.subTest(section=section), self.assertRaises(ValueError):
                resolved(source)


class NetCmcCompatibilityTests(unittest.TestCase):
    def test_build_gate_is_conditional_on_actual_activation(self):
        info = {'pass_max_version': 1, 'surface_adhesion_version': 1,
                'free_cmc_repulsion_version': 1, 'free_cmc_inner_repulsion_version': 1}
        cfg, _, _ = resolved()
        for version in (None, 0, 2):
            old = dict(info)
            if version is not None:
                old['cmc_net_potential_version'] = version
            with self.subTest(version=version), self.assertRaisesRegex(ValueError, 'net CMC.*Rebuild'):
                RUNNER.require_local_adhesion_build(cfg, old)
        RUNNER.require_local_adhesion_build(cfg, dict(info, cmc_net_potential_version=1))
        cfg['cmc']['free_g_L'] = START
        RUNNER.require_local_adhesion_build(cfg, info)
        cfg['cmc']['free_g_L'] = FULL
        cfg['cmc']['net_potential']['enabled'] = False
        RUNNER.require_local_adhesion_build(cfg, info)

    def test_restart_requires_complete_matching_net_fingerprint_and_version(self):
        cfg, meta, values = resolved()
        saved = checkpoint(cfg, meta, values)
        self.assertTrue(RUNNER.validate_restart(cfg, meta, saved, 1))
        for key in RUNNER.NET_CMC_CHECKPOINT_KEYS:
            broken = copy.deepcopy(saved)
            del broken['immutable_config'][key]
            with self.subTest(missing=key), self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, meta, broken, 1)
            broken = copy.deepcopy(saved)
            broken['immutable_config'][key] *= .5
            with self.subTest(changed=key), self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, meta, broken, 1)
        disabled = copy.deepcopy(cfg)
        disabled['cmc']['net_potential']['enabled'] = False
        old_cfg, old_meta, old_values = resolved(disabled)
        old_saved = checkpoint(old_cfg, old_meta, old_values)
        self.assertTrue(RUNNER.validate_restart(old_cfg, old_meta, old_saved, 1))
        with self.assertRaises(ValueError):
            RUNNER.validate_restart(cfg, meta, old_saved, 1)
        with self.assertRaises(ValueError):
            RUNNER.validate_restart(old_cfg, old_meta, saved, 1)
        # Even a zero-amplitude stray signature is not a legacy checkpoint.
        old_saved['immutable_config']['cmc_net_blend'] = 0
        with self.assertRaises(ValueError):
            RUNNER.validate_restart(old_cfg, old_meta, old_saved, 1)


@unittest.skipUnless(shutil.which('g++'), 'C++ parser integration requires g++')
class NetCmcCppWiringTests(unittest.TestCase):
    def test_derived_net_parameters_reach_cpp_for_both_endpoints_and_intermediate(self):
        with tempfile.TemporaryDirectory(prefix='net-cmc-config-') as temporary:
            directory = Path(temporary)
            source = directory/'parse.cpp'
            source.write_text('''#include "quickConfig.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
  const auto c=slurry::gr_re2::parseConfig(argc,argv);
  std::cout<<std::setprecision(17)<<c.cmc_net_blend<<" "
    <<c.cmc_net_contact_force<<" "<<c.cmc_net_barrier_force<<" "
    <<c.cmc_net_attraction_range<<" "<<c.cmc_net_repulsion_range<<" "
    <<c.cmc_net_reference_length<<" "<<c.free_cmc_repulsion_pressure<<" "
    <<c.free_cmc_inner_repulsion_work<<" "<<c.roughness_gap<<" "<<c.pass_max;
}
''')
            executable = directory/'parse'
            subprocess.run(['g++', '-std=c++17', '-I'+str(ROOT/'olb-1.9r0/src/slurry/gr_re2'),
                            str(source), '-o', str(executable)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for free in (START, (START+FULL)/2, FULL):
                config = case()
                config['cmc']['free_g_L'] = free
                _, meta, values = resolved(config)
                path = directory/'run.cfg'
                path.write_text(''.join('{}={}\n'.format(key, value) for key, value in values.items()))
                result = subprocess.run([str(executable), '--config', str(path)], check=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                parsed = list(map(float, result.stdout.split()))
                self.assertEqual(parsed[0], meta['cmc']['net_potential']['blend'])
                for i, key in enumerate(NET_KEYS[1:], start=1):
                    if key in values:
                        self.assertEqual(parsed[i], values[key])
                self.assertEqual(parsed[6:], [values['free_cmc_repulsion_pressure'],
                                              values['free_cmc_inner_repulsion_work'], 2e-9, 1])


if __name__ == '__main__':
    unittest.main()
