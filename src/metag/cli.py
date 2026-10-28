"""Command line: ``metag score`` (GPU) and the CPU tool commands.

    metag score [--only KEY] [--sampling-profile thorough] [--reactions FILE] [--out DIR]
    metag preflight INPUT -o REPORT
"""
import argparse
import os
import sys


def parse_args(argv=None):
    from metag.tools import command_help
    ap = argparse.ArgumentParser(prog="metag score", description="Score reactions (requires a GPU, UMA and xtb).",
                                 epilog=command_help(), formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reactions", default=os.environ.get("RXN_FILE"),
                    help="reactions JSON (default: RXN_FILE, else the packaged benchmark panel)")
    ap.add_argument("--only", action="append", default=None, help="score only this reaction key (repeatable)")
    ap.add_argument("--out", default=os.environ.get("METAG_OUT", "./metag_out"), help="output directory")
    ap.add_argument("--sampling-profile", choices=("current", "thorough"), default=None,
                    help="thorough uses stricter plateaus and larger pools; costs more, no accuracy guarantee")
    return ap.parse_args(argv)


def run_batch(reactions, keys, backend, out, log=None):
    """Score `keys`, checkpointing every outcome to OUT/metag_<tag>.json before the next reaction."""
    from metag.execution.filesystem import atomic_write_json
    from metag.pipeline import no_estimate, score_reaction
    log = log or (lambda s: print(s, flush=True))
    tag = keys[0] if len(keys) == 1 else "all"
    output_path = os.path.join(out, f"metag_{tag}.json")
    rows = []
    for key in keys:
        try:
            r = score_reaction(reactions[key], backend, key=key, log=log)
        except Exception as exc:
            r = no_estimate(key, reactions[key], f"{type(exc).__name__}: {exc}", log)
        rows.append(r)
        atomic_write_json(output_path, rows)
    log("\n==== MetaG ΔrG'° (kJ/mol, 298.15 K, pH 7, I = 0) ====")
    log(f"  {'reaction':14s} {'ΔG':>7s} {'U_samp':>7s} {'exp':>14s} {'err':>16s}")
    for r in rows:
        if r.get("dG") is None:
            log(f"  {r['reaction']:14s} NO ESTIMATE: {r.get('suspect', 'calculation failed')}")
            continue
        log(f"  {r['reaction']:14s} {r['dG']:7.1f} {r['U_samp']:7.1f} {str(r['exp']):>14s} {str(r['err']):>16s}")
    log(f"wrote {output_path}")
    return rows


def score_main(argv=None):
    a = parse_args(argv)
    from metag.reactions import load_reactions
    reactions = load_reactions(a.reactions)
    if not reactions:
        raise SystemExit(f"no reactions found in {a.reactions or 'the packaged panel'}")
    keys = a.only or list(reactions)
    unknown = [k for k in keys if k not in reactions]
    if unknown:
        raise SystemExit("unknown reaction key(s): " + ", ".join(unknown))
    from metag.settings import Settings
    settings = Settings.from_env(sampling_profile=a.sampling_profile)
    print(f"loading UMA {settings.model} ... reactions={len(keys)} profile={settings.sampling_profile}", flush=True)
    from metag.energetics.backend import UMABackend
    from metag.energetics.uma import load_uma
    backend = UMABackend(load_uma(settings.model), settings)
    os.makedirs(a.out, exist_ok=True)
    return run_batch(reactions, keys, backend, a.out)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    from metag.tools import COMMANDS, command_help, run_command
    if argv and argv[0] in COMMANDS:
        return run_command(argv[0], argv[1:])
    if argv and argv[0] in ("tools", "-h", "--help"):
        print(__doc__.strip() + "\n\n" + command_help())
        return
    if argv and argv[0] == "score":
        argv = argv[1:]
    return score_main(argv)


if __name__ == "__main__":
    main()
