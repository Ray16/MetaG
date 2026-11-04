import json
from pathlib import Path

from metag.tools.preflight import audit_reactions


def test_rare_sugar_candidates_have_structurally_admissible_declared_pools():
    path = Path(__file__).resolve().parents[2] / 'benchmark/candidates/itoh_1994_structures.json'
    reactions = json.loads(path.read_text())
    assert len(reactions) == 4
    assert all('exp' not in rx and set(rx['species_state_pools']) == set(rx['species'])
               for rx in reactions.values())
    report = audit_reactions(reactions)
    assert all(record['status'] == 'no_known_gap_detected' for record in report['reactions'].values())
