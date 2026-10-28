"""Runtime numerical settings, resolved once per scoring call.

Every setting that can change a species free energy lives here, and the species-cache key is built from
the same object that drives the computation, so a cached value can never describe a different
configuration. Environment variables are read when :meth:`Settings.from_env` is called, not at import.

The scientific model itself (routes, solvation model, protonation tables) is fixed; these settings only
control search effort, the electronic model checkpoint and the Hessian implementation.
"""
from __future__ import annotations

import os
from dataclasses import asdict, dataclass, replace

# Search-effort presets. "thorough" uses stricter plateaus and larger pools; it costs more and does not
# guarantee accuracy or complete sampling.
SAMPLING_PROFILES = {
    "current": {},
    "thorough": {"conv_tol": 0.5, "conv_hits": 3, "conv_max": 12, "sample_scale": 2.0, "keep_scale": 2.0},
}
SAMPLING_PROFILE_VERSION = "2026-10-11-search-effort-v1"

# Bump whenever species-level physics changes in a way the settings below do not capture (sampling,
# relaxation, thermal, solvation code). It is part of every species-cache key.
PHYSICS_VERSION = "2026-10-28-heavy-atom-dedup-v1"


@dataclass(frozen=True)
class Settings:
    model: str = "uma-s-1p2p1"          # registered UMA checkpoint
    conv_tol: float = 2.5               # kJ/mol: ensemble G and minimum E must move less than this ...
    conv_hits: int = 2                  # ... for this many consecutive seed batches
    conv_max: int = 8                   # cap on seed batches
    sample_scale: float = 1.0           # scales the per-batch conformer pool
    keep_scale: float = 1.0             # scales the number of conformers relaxed and solvated
    thermal_window_kj: float = 15.0     # provisional per-minimum RRHO window during sampling
    thermal_max: int = 10               # at most this many provisional RRHO corrections per batch
    thermal_validation_max: int = 64    # more sampled minima than this fails the species
    thermal_relax_steps: int = 400      # tight re-optimisation budget during minimum validation
    thermal_retry_steps: int = 1200     # one retry of an unresolved minimum with this budget
    seed_offset: int = 0                # conformer seed offset, for independent replicate searches
    hessian: str = "autograd"           # exact autograd UMA Hessian; "fd" = central differences
    dg_sanity_kj: float = 500.0         # |ΔG| above this x reaction extent is flagged (not rejected)
    sampling_profile: str = "current"

    def __post_init__(self):
        if self.hessian not in ("autograd", "fd"):
            raise ValueError(f"hessian must be 'autograd' or 'fd', got {self.hessian!r}")
        if self.sampling_profile not in SAMPLING_PROFILES:
            raise ValueError("sampling profile must be one of " + ", ".join(SAMPLING_PROFILES))
        if self.thermal_relax_steps < 1:
            raise ValueError("thermal_relax_steps must be positive")
        if self.seed_offset < 0:
            raise ValueError("seed_offset must be nonnegative")
        if self.conv_max < 1 or self.conv_hits < 1:
            raise ValueError("conv_max and conv_hits must be positive")

    @classmethod
    def from_env(cls, environ=None, *, sampling_profile=None) -> "Settings":
        """Resolve settings from the environment (``os.environ`` by default) at call time.

        A sampling profile is applied last and overrides the individual search-effort variables."""
        env = os.environ if environ is None else environ
        profile = (sampling_profile or env.get("METAG_SAMPLING_PROFILE", "current")).strip().lower()
        s = cls(
            model=env.get("UMA_MODEL", cls.model),
            conv_tol=float(env.get("CONV_TOL", cls.conv_tol)),
            conv_hits=int(env.get("CONV_HITS", cls.conv_hits)),
            conv_max=int(env.get("CONV_MAX", cls.conv_max)),
            sample_scale=float(env.get("SAMPLE_SCALE", cls.sample_scale)),
            keep_scale=float(env.get("KEEP_SCALE", cls.keep_scale)),
            thermal_window_kj=float(env.get("THERMAL_ENS_WIN", cls.thermal_window_kj)),
            thermal_max=int(env.get("THERMAL_ENS_N", cls.thermal_max)),
            thermal_validation_max=int(env.get("THERMAL_VALIDATION_MAX", cls.thermal_validation_max)),
            thermal_relax_steps=int(env.get("THERMAL_RELAX_STEPS", cls.thermal_relax_steps)),
            thermal_retry_steps=int(env.get("THERMAL_RETRY_STEPS", cls.thermal_retry_steps)),
            seed_offset=int(env.get("CONFORMER_SEED_OFFSET", cls.seed_offset)),
            hessian=env.get("UMA_HESSIAN", cls.hessian).strip().lower(),
            dg_sanity_kj=float(env.get("DG_SANITY_KJ", cls.dg_sanity_kj)),
        )
        return s.with_profile(profile)

    def with_profile(self, profile) -> "Settings":
        if profile not in SAMPLING_PROFILES:
            raise ValueError("sampling profile must be one of " + ", ".join(SAMPLING_PROFILES))
        return replace(self, sampling_profile=profile, **SAMPLING_PROFILES[profile])

    def as_dict(self) -> dict:
        return asdict(self)

    def species_key(self) -> dict:
        """Everything that changes a cached species free energy (the model is added by the backend)."""
        from metag.energetics.conformers import DEDUP_E_TOL
        key = {"solv": "alpb", "budget": "nrot-tiered-v1",
               "conv_tol": self.conv_tol, "conv_hits": self.conv_hits, "conv_max": self.conv_max,
               "sample_scale": self.sample_scale, "keep_scale": self.keep_scale,
               "physics": PHYSICS_VERSION, "hessian": self.hessian,
               "dedup": f"rmsd-v3-heavy-fallback-etol{DEDUP_E_TOL}",
               "thermal_ensemble": f"v3-w{self.thermal_window_kj:g}-n{self.thermal_max}",
               "thermal_validation_max": self.thermal_validation_max,
               "thermal_relax_steps": self.thermal_relax_steps,
               "thermal_retry_steps": self.thermal_retry_steps,
               "seed_offset": self.seed_offset}
        return key
