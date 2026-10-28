"""metag.settings.Settings: environment parsing, profile precedence, validation and the cache key."""
import dataclasses

import pytest

from metag.settings import SAMPLING_PROFILES, Settings


def test_from_env_parses_every_variable():
    env = {"UMA_MODEL": "uma-x", "CONV_TOL": "1.5", "CONV_HITS": "4", "CONV_MAX": "9",
           "SAMPLE_SCALE": "1.5", "KEEP_SCALE": "0.5", "THERMAL_ENS_WIN": "10", "THERMAL_ENS_N": "7",
           "THERMAL_VALIDATION_MAX": "32", "THERMAL_RELAX_STEPS": "200", "THERMAL_RETRY_STEPS": "900",
           "CONFORMER_SEED_OFFSET": "3", "UMA_HESSIAN": " FD ", "DG_SANITY_KJ": "400"}
    s = Settings.from_env(env)
    assert (s.model, s.conv_tol, s.conv_hits, s.conv_max) == ("uma-x", 1.5, 4, 9)
    assert (s.sample_scale, s.keep_scale, s.thermal_window_kj, s.thermal_max) == (1.5, 0.5, 10.0, 7)
    assert (s.thermal_validation_max, s.thermal_relax_steps, s.thermal_retry_steps) == (32, 200, 900)
    assert (s.seed_offset, s.hessian, s.dg_sanity_kj, s.sampling_profile) == (3, "fd", 400.0, "current")


def test_empty_environment_gives_defaults():
    assert Settings.from_env({}) == Settings()


def test_profile_overrides_individual_search_variables():
    env = {"CONV_TOL": "9", "CONV_HITS": "1", "SAMPLE_SCALE": "5", "METAG_SAMPLING_PROFILE": "current"}
    s = Settings.from_env(env, sampling_profile="thorough")
    assert s.sampling_profile == "thorough"
    assert (s.conv_tol, s.conv_hits, s.conv_max, s.sample_scale, s.keep_scale) == (0.5, 3, 12, 2.0, 2.0)
    # The env profile is used when no argument is given; "current" leaves the env values alone.
    assert Settings.from_env(dict(env, METAG_SAMPLING_PROFILE="Thorough")).conv_tol == 0.5
    assert Settings.from_env(env).conv_tol == 9.0
    assert Settings().with_profile("thorough").conv_max == 12


@pytest.mark.parametrize("kwargs", [{"hessian": "numeric"}, {"sampling_profile": "fast"},
                                    {"seed_offset": -1}, {"conv_max": 0}, {"thermal_relax_steps": 0}])
def test_invalid_settings_raise(kwargs):
    with pytest.raises(ValueError):
        Settings(**kwargs)


def test_invalid_env_hessian_and_profile_raise():
    with pytest.raises(ValueError):
        Settings.from_env({"UMA_HESSIAN": "numeric"})
    with pytest.raises(ValueError):
        Settings.from_env({}, sampling_profile="fast")
    with pytest.raises(ValueError):
        Settings().with_profile("fast")


def test_settings_are_frozen():
    with pytest.raises(dataclasses.FrozenInstanceError):
        Settings().conv_tol = 1.0


# Fields that never change a species free energy: the model is added to the key by the backend, the
# profile name is captured by the values it sets, and dg_sanity_kj is a reaction-level warning threshold.
_NOT_SPECIES_FIELDS = {"model", "sampling_profile", "dg_sanity_kj"}
_PERTURB = {"conv_tol": 1.0, "conv_hits": 5, "conv_max": 11, "sample_scale": 3.0, "keep_scale": 3.0,
            "thermal_window_kj": 7.0, "thermal_max": 3, "thermal_validation_max": 5,
            "thermal_relax_steps": 17, "thermal_retry_steps": 19, "seed_offset": 2, "hessian": "fd"}


def test_species_key_contains_every_result_affecting_field():
    fields = {f.name for f in dataclasses.fields(Settings)} - _NOT_SPECIES_FIELDS
    assert fields == set(_PERTURB)
    base = Settings().species_key()
    for name, value in _PERTURB.items():
        assert Settings(**{name: value}).species_key() != base, name
    assert base["solv"] == "alpb" and "physics" in base
    assert Settings(dg_sanity_kj=1.0).species_key() == base


def test_profile_changes_species_key_and_as_dict_roundtrips():
    s = Settings()
    assert s.with_profile("thorough").species_key() != s.species_key()
    assert Settings(**s.as_dict()) == s
    assert set(SAMPLING_PROFILES) == {"current", "thorough"}
