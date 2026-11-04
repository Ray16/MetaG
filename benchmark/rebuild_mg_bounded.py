"""Rebuild the archived Mg-bounded reference without reading MetaG predictions.

Requires the original source snapshot and eQuilibrator database; see benchmark/README.md.
"""
import argparse, hashlib, json, math, re, statistics
from collections import Counter, defaultdict
from pathlib import Path

from metag.tools import build_reference as R
from metag.validation import reference_transform as B
from metag.validation.reference_conditions import audit_conditions, number, magnesium_evidence, magnesium_total_activity_ambiguity

SALT = re.compile(r'\bc\s*\(\s*Mg(?:Cl2|SO4|\(OAc\)2)?\s*\)\s*=\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*(mM|M|mol)')


def mg_interval(row, meta):
    """(lo, hi) pMg_free bounds, or None if the row has no Mg issue (point pMg)."""
    pmg = number(row.get('p_mg'))
    if pmg is not None:
        return (pmg, 14.0) if magnesium_total_activity_ambiguity(row, meta) else None
    if not magnesium_evidence(row, meta):
        return (4.0, 14.0)
    for data in (row, meta or {}):
        for field in ('additional data', 'Cofactor', 'buffer'):
            m = SALT.search(str(data.get(field, '')))
            if m:
                c = float(m.group(1)) * (1e-3 if m.group(2) == 'mM' else 1)
                if c > 0:
                    return (max(-math.log10(c), 0.0), 14.0)
    return (2.0, 14.0)


def base_exclusion(row, meta):
    """The 10-10 generic exclusions, minus the three Mg rules handled by bounding."""
    why = R.measurement_exclusion(row, meta)
    if why in ('Mg_condition_unspecified', 'Mg_present_free_activity_unknown',
               'pMg_matches_total_magnesium_recipe'):
        # The generic checker returns early for these Mg conditions; after
        # replacing its Mg exclusion with bounding, retain its later I check.
        ionic = number(row.get('ionic_strength'))
        return 'invalid_ionic_strength' if ionic is not None and ionic < 0 else None
    return why


def same_chemistry(hits, ms):
    def full(rid, o):
        return tuple(sorted((c, round(o * v, 6)) for c, v in ms[rid].items() if c != B.PROTON))
    return len({full(r, o) for r, o in hits}) == 1


def rebuild(inputs, threshold, tr, ms, aliases, index, rows, *, excluded_records=None):
    used, counts = defaultdict(list), Counter()
    for row, equation, how in rows:
        st, ok = B.parse_equation(equation)
        if not ok or any(c not in aliases for c in st):
            continue
        mapped = defaultdict(float)
        for c, k in st.items():
            mapped[aliases[c]] += k
        hits = sorted(set(index.get(B.signature(mapped, {B.PROTON, B.WATER}), [])))
        panel = [h for h in hits if h[0] in inputs]
        if not panel:
            continue
        rid, orient = panel[0]
        raw = {k: v for k, v in row.items() if not k.startswith('_')}
        meta = row['_reference_metadata']
        why = base_exclusion(raw, meta)
        if why is None and len(hits) != 1 and not same_chemistry(hits, ms):
            why = 'ambiguous_reaction_mapping'
        if why is None and tr.reaction(rid) is None:
            why = 'missing_balanced_transformable_reaction'
        rec = None
        if why is None:
            t, ph = number(raw['temperature']), number(raw['p_h'])
            ionic = number(raw.get('ionic_strength'))
            ionic = 0.25 if ionic is None else ionic
            observed = orient * -B.R_KJ * t * math.log(number(raw['K_prime']))
            iv = mg_interval(raw, meta)
            grid = [number(raw['p_mg'])] if iv is None else [iv[0] + (iv[1] - iv[0]) * i / 10 for i in range(11)]
            try:
                vals = [observed + tr.correction(rid, ph, ionic, g, t) for g in grid]
            except (ValueError, RuntimeError, KeyError, Warning) as exc:
                vals, why = None, 'condition_transform_failed'
            if vals:
                spread = max(vals) - min(vals)
                if spread > threshold:
                    why = 'Mg_sensitive_beyond_threshold'
                else:
                    rec = dict(reaction=rid, source_row=row['_source_row'], T=t, pH=ph, I=ionic,
                               pMg_interval=iv, mg_spread=spread, reference_dG=(max(vals) + min(vals)) / 2,
                               ionic_imputed=number(raw.get('ionic_strength')) is None)
        if why:
            counts[why] += 1
            if excluded_records is not None:
                excluded_records.append(dict(reaction=rid, source_row=row['_source_row'],
                                             equation=equation, mapping=how,
                                             orientation=orient, mapping_candidates=hits,
                                             reason=why))
        else:
            used[rid].append(rec)
    labels = {rid: dict(exp=round(statistics.median(r['reference_dG'] for r in recs), 2), n=len(recs),
                        sd=round(statistics.stdev([r['reference_dG'] for r in recs]), 2) if len(recs) > 1 else None,
                        max_mg_spread=round(max(r['mg_spread'] for r in recs), 3),
                        temperatures=sorted({r['T'] for r in recs}),
                        n_298=sum(abs(r['T'] - 298.15) < 1e-6 for r in recs),
                        I_imputed=any(r['ionic_imputed'] for r in recs))
              for rid, recs in used.items()}
    return labels, dict(counts), used


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reactions', type=Path, required=True, help='original frozen prediction_structures.json')
    parser.add_argument('--source', type=Path, required=True, help='original openTECR CSV snapshot directory')
    parser.add_argument('--modelseed', type=Path, required=True, help='original ModelSEED snapshot directory')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True, help='measurement-level records and input provenance')
    args = parser.parse_args(argv)
    inputs = json.loads(args.reactions.read_text())
    ms = B.load_modelseed_reactions(args.modelseed)
    aliases = B.kegg_to_modelseed(args.modelseed)
    index = defaultdict(list)
    for rid, st in ms.items():
        index[B.signature(st, {B.PROTON, B.WATER})].append((rid, 1))
        index[B.signature({c: -v for c, v in st.items()}, {B.PROTON, B.WATER})].append((rid, -1))
    rows, mapping = R.opentecr_rows(args.source)
    tr = B.Transformer(ms)
    out, measurements, exclusions = {}, {}, {}
    for thr in (2.0, 1.0, 4.0):
        excluded_rows = []
        labels, counts, used = rebuild(inputs, thr, tr, ms, aliases, index, rows,
                                       excluded_records=excluded_rows)
        out[str(thr)] = dict(n_labels=len(labels), n_measurements=sum(len(v) for v in used.values()),
                             excluded=counts, labels=labels)
        measurements[str(thr)] = dict(used)
        exclusions[str(thr)] = excluded_rows
        print(thr, len(labels), out[str(thr)]['n_measurements'], counts, flush=True)
    paths = [args.reactions, *sorted(args.source.glob('*.csv')),
             args.modelseed / 'Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt',
             *sorted((args.modelseed / 'Biochemistry').glob('reaction_??.tsv'))]
    provenance = {'reaction_sha256': hashlib.sha256(args.reactions.read_bytes()).hexdigest(),
                  'source_sha256': {str(p.relative_to(args.source)): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in paths if p.parent == args.source},
                  'modelseed_sha256': {str(p.relative_to(args.modelseed)): hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in paths if args.modelseed in p.parents},
                  'transform': tr.provenance(), 'mapping_summary': dict(mapping),
                  'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'note': 'Reproduction of the historical policy, including its original assumptions.'}
    for path, data in ((args.output, out), (args.audit, {'provenance': provenance,
                                                        'used_measurements': measurements,
                                                        'excluded_panel_measurements': exclusions})):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
