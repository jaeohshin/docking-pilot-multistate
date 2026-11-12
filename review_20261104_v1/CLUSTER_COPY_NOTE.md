# Cluster copy of this review — omitted evidence archives

This review was produced on 2026-11-04 on the desktop machine (`/data/work/sidechain-multistate`,
parent HEAD `193c6c5`), reading the completed docking study on the cluster over SSH. It was committed
to the cluster checkout on 2026-11-12, because the cluster is where this study's data and history live.

Two evidence archives in `evidence/` are **not** copied here, because on the cluster they would duplicate
their own sources in the same tree:

| Omitted file | Bytes | SHA256 (first 16) | Source on the cluster |
|---|---:|---|---|
| `evidence/docking_source_8b6814b.tar.gz` | 9,520,509 | `c4625035d90372d6` | this repository at commit `8b6814b` (`git archive 8b6814b`) |
| `evidence/selected_existing.tar.gz` | 13,196,082 | `a530d59954eb8ca3` | the files listed in `evidence/COLLECTED_MANIFEST.json`, each with its own SHA256 |

Both full-archive hashes are also recorded in `SHA256SUMS.json`, whose `total_bytes` (34,998,120) refers to
the complete desktop folder including them. The committed tree here is 12,304,651 bytes in 134 files.
Nothing else was changed, added or regenerated: file contents are byte-identical to the desktop originals.

The desktop copy, including both archives, is preserved at
`/data/work/sidechain-multistate/experiments/docking_pilot_20261028/review_20261104_v1/`.
It is untracked there. If that machine is lost, `docking_source_8b6814b.tar.gz` is reproducible from this
repository's history, and `selected_existing.tar.gz` from the paths and hashes in `COLLECTED_MANIFEST.json`
as long as the `expand_v1`/`followup_v1` outputs remain.
