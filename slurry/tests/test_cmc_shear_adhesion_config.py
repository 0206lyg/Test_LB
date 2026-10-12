"""Unified CMC concentration law, steady-shear scaling and restart safety."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'shear_cmc_runner', ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
START = 1.1015785714285715
FULL = 13.101578571428572


def source():
    return json.loads((ROOT/'slurry/cases/gr_CMC.json').read_text())


def resolved(config=None, **overrides):
    cfg, meta = RUNNER.resolve(source() if config is None else config, **overrides)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    return cfg, meta, values


def checkpoint(cfg, meta, values):
    immutable = {}
    for physical, version, keys in (
            ('surface_adhesion', 'interaction_model_version', RUNNER.SURFACE_CHECKPOINT_KEYS),
            ('cmc_net_blend', 'cmc_net_potential_version', RUNNER.NET_CMC_CHECKPOINT_KEYS),
            ('current_adhesion_rolling', 'current_adhesion_rolling_version',
             RUNNER.CURRENT_ADHESION_ROLLING_CHECKPOINT_KEYS)):
        if values.get(physical, 0) > 0:
            immutable.update({key: values[key] for key in keys if key != version})
            immutable[version] = 1
    immutable.update(dt_s=meta['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    return {'step': 1, 'immutable_config': immutable}


class ShearCmcConstitutiveTests(unittest.TestCase):
    def test_shipped_case_has_one_concentration_law_and_unchanged_contact_parameters(self):
        config = source()
        self.assertNotIn('net_potential', config['cmc'])
        self.assertNotIn('free_repulsion', config['cmc'])
        self.assertTrue(config['cmc']['shear_adhesion']['enabled'])
        self.assertEqual(config['cmc']['free_g_L'], FULL)
        cfg, meta, values = resolved(config)
        self.assertEqual(values['current_adhesion_rolling'], 1)
        self.assertEqual(values['sliding_friction'], .1)
        self.assertEqual(values['tangential_stiffness'], 80)
        self.assertEqual(values['rolling_length'], 1e-7)
        self.assertEqual(values['rolling_yield_angle'], .01)
        self.assertEqual(cfg['fluid']['dynamic_viscosity_Pa_s'], .00089)
        self.assertTrue(meta['cmc']['shear_adhesion']['active'])

    def test_concentration_endpoints_and_saturation_never_restore_legacy_repulsion(self):
        for free, expected_p, expected_retention in (
                (0, 0, 1), (START, 0, 1), ((START+FULL)/2, .5, .55),
                (FULL, 1, .1), (FULL+12, 1, .1)):
            config = source()
            config['cmc']['free_g_L'] = free
            _, meta, values = resolved(config, shear_rate=10)
            shear = meta['cmc']['shear_adhesion']
            with self.subTest(free=free):
                self.assertEqual(shear['p'], expected_p)
                self.assertAlmostEqual(shear['retention'], expected_retention)
                self.assertEqual(values['cmc_net_blend'], 1)
                self.assertEqual(values['cmc_net_barrier_force'], 0)
                self.assertEqual(meta['cmc']['net_potential']['re2_tail_weight'], 0)
                self.assertFalse(any(key.startswith('free_cmc_') for key in values))

    def test_reference_force_is_derived_from_adsorption_and_not_independent_force_input(self):
        cfg, meta, values = resolved()
        interaction = cfg['interaction']
        h0 = cfg['rough_contact']['roughness_gap_m']
        bare_work = interaction['adhesion_work_J_m2']
        background = interaction['hamaker_J']/(12*math.pi*h0*h0) * (
            1-(interaction['sigma_lj_m']/h0)**6/30)
        expected_work = background + .33*(bare_work-background)
        a, c = cfg['particles']['diameter_m']/2, cfg['particles']['thickness_m']/2
        pure_reference_force = math.pi*a*a/c*bare_work
        expected_force = math.pi*a*a/c*expected_work
        shear = meta['cmc']['shear_adhesion']
        self.assertAlmostEqual(shear['baseline_contact_force_N']/expected_force, 1)
        self.assertAlmostEqual(values['cmc_net_contact_force']/expected_force, 1)
        self.assertAlmostEqual(expected_force/pure_reference_force, .3501, places=4)
        self.assertGreater(expected_force, 3.27e-7)
        self.assertLess(expected_force, 3.29e-7)

    def test_rate_override_changes_only_effective_strength_and_recovers_high_rate_limit(self):
        config = source()
        original = copy.deepcopy(config)
        baseline = resolved(config)[1]['cmc']['shear_adhesion']['baseline_contact_force_N']
        for rate, retention in ((1, .01), (10, .1), (100, 1), (200, 1)):
            cfg, meta, values = resolved(config, shear_rate=rate)
            shear = meta['cmc']['shear_adhesion']
            with self.subTest(rate=rate):
                self.assertEqual(cfg['flow']['shear_rate_s_inv'], rate)
                self.assertAlmostEqual(shear['retention'], retention)
                self.assertEqual(shear['baseline_contact_force_N'], baseline)
                self.assertAlmostEqual(values['cmc_net_contact_force']/baseline, retention)
                self.assertEqual(shear['effective_contact_force_N'], values['cmc_net_contact_force'])
                self.assertEqual(values['cmc_net_attraction_range'], 1e-9)
                self.assertEqual(values['cmc_net_barrier_force'], 0)
        self.assertEqual(config, original)

    def test_low_free_concentration_retains_the_same_baseline_at_every_rate(self):
        config = source()
        config['cmc']['free_g_L'] = START
        forces = []
        for rate in (1, 10, 100):
            _, meta, values = resolved(config, shear_rate=rate)
            self.assertEqual(meta['cmc']['shear_adhesion']['retention'], 1)
            forces.append(values['cmc_net_contact_force'])
        self.assertEqual(forces, [forces[0]]*3)
        self.assertEqual(forces[0], resolved()[2]['cmc_net_contact_force'])

    def test_adsorption_clamp_does_not_transfer_mass_to_free_cmc(self):
        config = source()
        config['cmc']['free_g_L'] = (START+FULL)/2
        _, meta, baseline_values = resolved(config, shear_rate=10)
        config['cmc']['adsorbed_g_L'] *= 3
        _, clamped, values = resolved(config, shear_rate=10)
        self.assertTrue(clamped['cmc']['adsorbed_clamped'])
        self.assertEqual(clamped['cmc']['free_g_L'], (START+FULL)/2)
        self.assertEqual(clamped['cmc']['shear_adhesion'], meta['cmc']['shear_adhesion'])
        self.assertEqual(values['cmc_net_contact_force'], baseline_values['cmc_net_contact_force'])

    def test_resolving_saved_config_does_not_scale_twice_and_geometry_is_recomputed(self):
        cfg, meta, values = resolved(shear_rate=10)
        roundtrip = resolved(json.loads(json.dumps(cfg)))
        self.assertEqual(roundtrip, (cfg, meta, values))
        cfg['particles']['thickness_m'] *= 2
        new_values = resolved(cfg)[2]
        self.assertAlmostEqual(new_values['cmc_net_contact_force']/values['cmc_net_contact_force'], .5)
        self.assertAlmostEqual(new_values['cmc_net_reference_length']/values['cmc_net_reference_length'], .5)

    def test_zero_shear_constitutive_limit_is_separate_from_positive_flow_validation(self):
        cfg, _, _ = resolved()
        cfg['flow']['shear_rate_s_inv'] = 0
        _, cmc_meta = RUNNER.resolve_cmc(cfg)
        self.assertEqual(cmc_meta['shear_adhesion']['retention'], 0)
        self.assertEqual(cmc_meta['shear_adhesion']['effective_contact_force_N'], 0)
        with self.assertRaises(ValueError):
            RUNNER.resolve(cfg)


class ShearCmcValidationTests(unittest.TestCase):
    def test_invalid_shear_law_is_rejected_before_running(self):
        for bad in (None, [], True, 'enabled'):
            config = source()
            config['cmc']['shear_adhesion'] = bad
            with self.subTest(section=bad), self.assertRaises(ValueError):
                resolved(config)
        for updates in ({'enabled': 1}, {'unknown': 1}, {'start_free_g_L': -1},
                        {'full_free_g_L': START}, {'full_shear_rate_s_inv': 0},
                        {'full_shear_rate_s_inv': math.inf}, {'attraction_range_m': 0},
                        {'attraction_range_m': math.nan}):
            config = source()
            config['cmc']['shear_adhesion'].update(updates)
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                resolved(config)

    def test_enabled_legacy_models_cannot_be_silently_combined_with_shear_law(self):
        for key in ('net_potential', 'free_repulsion'):
            config = source()
            config['cmc'][key] = {'enabled': True}
            with self.subTest(section=key), self.assertRaises(ValueError):
                resolved(config)

    def test_new_mode_requires_contact_and_current_adhesion_rolling(self):
        for section, key in (('rough_contact', 'enabled'),
                             ('rough_contact', 'current_adhesion_rolling'),
                             ('interaction', 'surface_adhesion')):
            config = source()
            config[section][key] = False
            with self.subTest(section=section, key=key), self.assertRaises(ValueError):
                resolved(config)

    def test_build_capability_gate_applies_to_new_mode_only(self):
        cfg, _, _ = resolved()
        info = {'pass_max_version': 1, 'surface_adhesion_version': 1,
                'cmc_net_potential_version': 1, 'current_adhesion_rolling_version': 1}
        for version in (None, 0, 2):
            old = dict(info)
            if version is not None:
                old['cmc_shear_adhesion_version'] = version
            with self.subTest(version=version), self.assertRaisesRegex(ValueError, 'shear.*Rebuild'):
                RUNNER.require_local_adhesion_build(cfg, old)
        RUNNER.require_local_adhesion_build(cfg, dict(info, cmc_shear_adhesion_version=1))
        pure = json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())
        RUNNER.require_local_adhesion_build(RUNNER.resolve(pure)[0], info)

    def test_restart_rejects_changed_adhesion_amplitude(self):
        cfg, meta, values = resolved(shear_rate=10)
        saved = checkpoint(cfg, meta, values)
        self.assertTrue(RUNNER.validate_restart(cfg, meta, saved, 1))
        changed = copy.deepcopy(cfg)
        changed['cmc']['free_g_L'] = START
        changed_cfg, changed_meta, _ = resolved(changed)
        with self.assertRaisesRegex(ValueError, 'cmc_net_contact_force'):
            RUNNER.validate_restart(changed_cfg, changed_meta, saved, 1)
        broken = copy.deepcopy(saved)
        broken['immutable_config']['cmc_net_barrier_force'] = 6e-11
        with self.assertRaises(ValueError):
            RUNNER.validate_restart(cfg, meta, broken, 1)


if __name__ == '__main__':
    unittest.main()
