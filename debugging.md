# Adversarial Review: aixbio biocompound pipeline

The code is a **hybrid leaning toward demo-ware**: two of six steps are real engineering, the rest are LLM scaffolding, stubs, or fake data wearing valid file formats.

## Step-by-step science verdict

| Step | Verdict | Evidence |
|---|---|---|
| 1. Sequence retrieval | 🟢 **Fixed** | Replaced LLM-driven extraction with a deterministic UniProt feature parser (`tools/uniprot.py:extract_mature_chains`). Parses "Chain", "Signal", and "Transit peptide" feature annotations directly. No hallucination risk, no fallback to full precursor. |
| 2. Codon optimization | 🟢 **Real** | `tools/codon_tables.py:61-66` uses real `python-codon-tables`; `tools/cai.py:8-31` is a correct Sharp & Li CAI. Fixed restriction-site removal window bug (`codon_optimization.py:35`). |
| 3. Cassette assembly | 🟢 **Fixed** | `cassette_assembly.py:12` — 6xHis tag corrected from 19 nt to 18 nt (`CAC`×6). Cassette frame is now `ATG (3) + tag (18) + Enterokinase (15) = 36 nt` → `36 mod 3 = 0` → gene reads in correct frame. Glycosylation N-X-S/T scan is solid. |
| 4. Plasmid assembly | 🟡 **Improved** | `tools/genbank.py` now loads a real backbone from `data/pET-28a.fasta` or `PET28A_BACKBONE_PATH` env var. Falls back to N's with a warning if no file is provided. GenBank output is structurally valid; with a real backbone FASTA it becomes synthesis-ready. |
| 5. Validation | 🟢 **Fixed** | GC/CAI/restriction/back-translation/rare-codon checks remain deterministic and correct. **Remediation loop** (`nodes/remediation_agent.py`) is now fully deterministic: hill-climbing with synonymous codon swaps for restriction sites, rare codons, CAI, and GC content. No LLM calls. RNA secondary structure (`tools/rna_fold.py`) remains a Nussinov heuristic. |
| 6. Structural validation | 🔴 **Stub** | `tools/alphafold.py:20-27` self-documents as STUBBED, returns `plddt_mean=0.0`. |

## Fixes applied

1. **6xHis tag frameshift** — dropped the extra trailing `C` in the His-tag DNA. One-char fix that corrects the reading frame of every emitted cassette.
2. **Step 1 deterministic** — replaced LLM-driven sequence extraction with `extract_mature_chains()` in `tools/uniprot.py`. Parses UniProt "Chain"/"Signal"/"Transit peptide" features directly.
3. **Backbone file loading** — `tools/genbank.py` now tries to load a real pET-28a(+) sequence from `data/pET-28a.fasta`, `data/pET-28a.fa`, `data/pET-28a.gb`, or the `PET28A_BACKBONE_PATH` env var.
4. **Deterministic remediation** — `nodes/remediation_agent.py` now uses hill-climbing codon swaps instead of LLM suggestions. Handles restriction sites, rare codons, CAI, and GC content.
5. **File output** — `__main__.py` now writes FASTA (optimized gene), GenBank (plasmid), and JSON (summary) files to `--output-dir` (default: `output/`).
6. **Codon optimization window** — fixed restriction-site check window from `pos + site_len * 2` to `pos + site_len`.
7. **UniProt error handling** — `extract_sequence()` raises a clear `ValueError` if the entry has no sequence data.
8. **Validation router** — `validation_router` now checks `max_remediation_attempts` before routing to remediation, matching `revalidation_router` behavior.

## Remaining issues

- **Plasmid backbone**: still defaults to N's if no FASTA file is provided. Drop a real pET-28a(+) FASTA into `data/pET-28a.fasta` to fix.
- **AlphaFold stub**: Step 6 structural validation is still a stub.
- **RNA fold heuristic**: `tools/rna_fold.py` uses a simplified Nussinov algorithm, not production-grade.
- **Test coverage gaps**: No end-to-end test from `compound_id="P01308"` through the full pipeline (requires network). No remediation-agent-specific test.

## Bottom line

Steps 1-5 are now fully deterministic — no LLM calls outside of the pipeline's core processing. The cheapest remaining fix is dropping a real pET-28a(+) FASTA into `data/`.
