# Building a publication review snapshot

`tools/publication_snapshot.py` packages two explicitly named Git commits: the
current MetaG code and the separate manuscript repository. It reads committed
trees, so later commits, uncommitted edits and untracked files cannot silently
enter the selected snapshot. It preserves file permissions and records each
member's SHA-256 and source. An existing output directory is never overwritten.

The legacy code/evidence archive supplies historical raw studies and notices.
Its old runtime is retained under `historical_snapshots/review_20261104/`, away
from the current importable `MetaG/src` tree. Old distribution and verification
records remain identified as historical. The three companion archives retain
identical bytes and their original study layouts. No historical calculation is
assigned the identity of current development code.

From a MetaG checkout:

```sh
python tools/publication_snapshot.py \
  --code-repo /path/to/MetaG \
  --paper-repo /path/to/MetaG/manuscript \
  --code-commit FULL_CODE_COMMIT \
  --paper-commit FULL_MANUSCRIPT_COMMIT \
  --legacy-dir /path/to/verified/previous/review \
  --output /path/to/new/review
```

The legacy directory must contain `UPLOAD_MANIFEST.json` and the four named
code/study/comparison archives. Their sizes and hashes are checked before
packaging. New archives receive complete member readback checks. Already
compressed evidence is stored without recompression. The builder uses one CPU
process and never launches molecular calculations or retrieves model weights.

`SNAPSHOT_MANIFEST.json` records both commits, every new archive member, the
unchanged companion hashes and the snapshot's limited scope. The source archive
and manuscript archive contain no Git internals. Extract all four code/study
archives into the same parent to retain the historical study layout.

The extracted replay can be run with the supplied Linux verifier:

```sh
python tools/verify_publication_snapshot.py \
  --snapshot /path/to/new/review \
  --destination /path/to/new/extraction \
  --logs /path/to/new/check-logs \
  --forbid-root /path/to/original/project
```

The verifier authenticates the input archive hashes, extracts regular files,
runs every manuscript checker in a fresh Python process and verifies source
hashes again afterwards. An inherited Python read guard blocks reads under the
original project roots and network connections. Guard receipts are required for
each direct checker; even a caught fallback attempt fails the portable result.
The supplied tests exercise that failure case as well as path-escape rejection.
Native system calls are not an operating-system sandbox. Use the declared CPU
audit environment; the verifier does not install dependencies or run quantum jobs.

Compiled PDFs and build/check receipts require a separate verification against
the selected manuscript commit. Run the manuscript checks from the extracted
snapshot with `METAG_EVIDENCE_PORTABLE_ONLY=1`; access to original project files
must be excluded before claiming portable verification. Do not infer scientific
accuracy from archive integrity or CPU checks.

The current package is a review snapshot. Independent experimental validation,
selection of the final supported production reference policy, representative
numerical ModelSEED validation and user-deferred submission metadata remain
unfinished. A public archival deposit and final distribution identity also need
to correspond to the eventual validated method; this builder does not publish
a release or submit a manuscript.
