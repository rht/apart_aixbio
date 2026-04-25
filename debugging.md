# Adversarial Review: AixBio Pipeline

## Verdict

Competent scaffolding, but several bugs and fundamental shortcuts that would produce incorrect biological output.

---

## CRITICAL: Bugs That Produce Wrong Results

### 1. 6xHis tag encodes the wrong protein

**File:** `aixbio/nodes/cassette_assembly.py:7`

```python
TAG_SEQUENCES = {"6xHis": "CACCACCACCACCACCAC"}
```

This is 17 nucleotides -- not divisible by 3. It will **frameshift the entire downstream gene**. A 6xHis tag is 6 histidines = `CAC CAC CAC CAC CAC CAC` = 18 nucleotides. The sequence here is missing a `C`. This single typo makes every cassette produced by the pipeline encode a completely wrong protein.

### 2. Validation runs on the wrong DNA

**File:** `aixbio/nodes/sequence_validation.py:16`

```python
dna = dna_chain.dna_sequence  # the raw optimized gene
```

But the actual construct that goes into the host cell is the **cassette** (`ATG + tag + protease + gene + stop`). The validation checks GC content, CAI, restriction sites, and RNA secondary structure on just the gene -- not the actual expression construct. Restriction sites in the tag or protease sequences would go undetected. The 5' RNA structure check is particularly wrong since the ribosome encounters the tag region first, not the gene.

### 3. Codon optimization restriction-site removal is fragile

**File:** `aixbio/nodes/codon_optimization.py:26-29`

```python
for ci in range(codon_start, min(codon_end, len(codons))):
    alts = synonymous_alternatives(codons[ci])
    if alts:
        codons[ci] = alts[0]
        break  # only swaps ONE codon per site
```

It picks `alts[0]` -- the first alphabetical alternative, not one that actually eliminates the restriction site. If the first alternative still contains part of the recognition sequence (which is likely since restriction sites span codon boundaries), the site persists. The 5-pass loop provides some resilience, but there's no guarantee the greedy choice converges.

### 4. Enterokinase recognition site is wrong

**File:** `aixbio/nodes/cassette_assembly.py:11`

```python
PROTEASE_SITE_SEQUENCES = {"Enterokinase": "GATGATGATGATAAAG"}
```

Enterokinase cleaves after DDDDK (Asp-Asp-Asp-Asp-Lys). The DNA for DDDDK would be `GAT GAT GAT GAT AAA` or `AAG` for the lysine. The sequence here is `GATGATGATGATAAAG` = 16 nt, which is not divisible by 3 either. Combined with the 6xHis bug, the frame is now doubly shifted.

### 5. CNBr "protease site" is just ATG

**File:** `aixbio/nodes/cassette_assembly.py:13`

```python
"CNBr": "ATG"
```

CNBr is a chemical cleavage reagent that cuts after methionine, not a protease with a recognition sequence you encode. Using `ATG` (Met) as the site means the cassette would just have `ATG ATG` (start + "site") which adds an extra methionine. This is misleading at best.

---

## MAJOR: Architectural & Correctness Issues

### 6. `translate_dna` silently drops trailing nucleotides

**File:** `aixbio/tools/codon_tables.py:52`

```python
return "".join(translate_codon(c) for c in codons if len(c) == 3)
```

If the DNA length isn't divisible by 3 (e.g., after the frameshift from bug #1), the trailing 1-2 nucleotides are silently dropped. This masks frameshifts rather than catching them.

### 7. RNA fold estimator is biologically meaningless

**File:** `aixbio/tools/rna_fold.py:18-23`

The heuristic pairs the first nucleotide with the last, second with second-to-last, etc. -- a simple palindrome matcher. Real RNA secondary structure involves non-sequential base pairing, loops, bulges, and stacking energies. This will false-negative on actual stable hairpins and false-positive on palindromic sequences. The -10 kcal/mol threshold is calibrated for real folding algorithms; applying it to this heuristic produces meaningless results.

### 8. GenBank plasmid is fake

**File:** `aixbio/tools/genbank.py:22`

```python
plasmid_seq = "N" * PET28A_BACKBONE_SIZE + cassette_dna
```

The "backbone" is just 5369 N's. No promoter, no RBS, no terminator, no origin of replication, no antibiotic resistance gene. The GenBank file is structurally valid but biologically useless -- you couldn't order this and have it work. The `pET28A_BACKBONE_SIZE = 5369` is stated as fact but never verified.

### 9. Plasmid assembly adds restriction sites that Step 2 removed

**File:** `aixbio/nodes/plasmid_assembly.py:14-16`

```python
flank_5 = ENZYME_SITES.get(cloning_sites[0], "")  # e.g., GGATCC for BamHI
flank_3 = ENZYME_SITES.get(cloning_sites[1], "")  # e.g., CTCGAG for XhoI
flanked_dna = flank_5 + cassette.full_dna + flank_3
```

Step 2 removes BamHI and XhoI sites from the gene. Step 4 adds them right back at the flanks. If validation were run on the flanked cassette (which it should be per point #2), it would always fail the restriction sites check. The current code avoids this by validating the un-flanked gene -- but that's the wrong thing to validate.

### 10. AlphaFold is entirely stubbed

**File:** `aixbio/tools/alphafold.py`

Returns all zeros. Step 6 will always report pLDDT=0, RMSD=0, perplexity=0. Any downstream use of these signals is worthless. This isn't acknowledged prominently in the CLI output.

---

## MODERATE: Design Concerns

### 11. Spec says "Composable Step Chain" but implementation uses LangGraph

The `handoff.md` describes a simple functional pipeline (`pipeline()`, `retry()`, `fork()` -- ~50 LOC). The actual implementation uses LangGraph with `StateGraph`, `Send`, `interrupt()`, `MemorySaver`, etc. -- a heavy framework dependency. This isn't inherently wrong, but it's a complete departure from the spec's design rationale section. The stated advantages of simplicity and debuggability are lost.

### 12. LLM is in the critical path for chain extraction with a fragile fallback

**File:** `aixbio/nodes/sequence_retrieval.py:48`

If the LLM response can't be parsed, the fallback uses the **full precursor sequence** including signal peptides and pro-domains:

```python
def _fallback_single_chain(compound_id, name, sequence):
    return {"chains": [{"id": f"{name}_{compound_id}", "aa_sequence": sequence, ...}]}
```

This is the exact opposite of what you want. For insulin, the fallback would produce a 110-aa preproinsulin instead of the 21-aa A-chain and 30-aa B-chain. The pipeline would happily optimize and assemble a protein that cannot fold correctly.

### 13. Glycosylation warning detection is brittle

**File:** `aixbio/nodes/cassette_assembly.py:42-44`

```python
chain_id_upper = dna_chain.id.upper()
for compound in glycosylation_compounds:
    if compound.upper() in chain_id_upper:
```

This checks if "EPO" or "TPA" appears in the chain ID string. If the LLM names the chain "Erythropoietin_A" or "tPA_Kringle", it won't match "EPO" or will coincidentally match "TPA". This should be based on the protein record / UniProt annotations, not string matching on an LLM-generated ID.

### 14. No retry/timeout on the UniProt API call

`fetch_uniprot_entry` makes a single HTTP call with a 30s timeout. No retries. UniProt has intermittent availability issues. A single network hiccup kills the pipeline.

### 15. `asyncio.run()` called inside a sync node

**File:** `aixbio/nodes/sequence_retrieval.py:21`

```python
entry = asyncio.run(fetch_uniprot_entry(compound_id))
```

If LangGraph is already running an event loop (e.g., async execution), this will raise `RuntimeError: cannot run nested event loop`. This works now because the graph is invoked synchronously, but it's a landmine for any future async execution context.

### 16. Rare codon list includes CUA but the codon table uses uppercase T not U

**File:** `aixbio/tools/codon_tables.py:43`

`RARE_CODONS_ECOLI` contains `"CUA"` but DNA codons use `T` not `U`. The check at `sequence_validation.py:63` compares `c.upper()` against the set -- but `CUA` is an RNA codon, never produced by the DNA-based codon optimizer. This rare codon will never be detected. Should be `CTA`.

---

## MINOR: Test Coverage & Code Quality

### 17. Tests only cover the happy path for one compound

All tests use insulin B-chain (30 aa). No tests for multi-chain proteins, long sequences, edge-case amino acids (selenocysteine), or the remediation loop. The LLM-dependent nodes (Step 1, Step 5b) have zero test coverage.

### 18. `remediation_history` uses `append_log` reducer but `apply_fixes` returns a plain list

In `chain_state.py:36`, `remediation_history` uses `Annotated[list[RemediationAction], append_log]`, so each return appends. But `apply_fixes` returns `{"remediation_history": actions}` where `actions` is the list from *this round only*. This means the history accumulates correctly via the reducer, but it's non-obvious and easy to break.

### 19. `package_result_failed` is identical to `package_result`

**File:** `aixbio/nodes/merge_results.py:24`

```python
def package_result_failed(state):
    return package_result(state)
```

Failed chains are packaged identically to passing chains. There's no status field that distinguishes "passed validation" from "gave up after max retries." The `validation_passed` field reflects the last check, which will be `False`, but there's no "exceeded max attempts" flag.

---

## Summary Scorecard

| Area | Rating | Notes |
|------|--------|-------|
| Architecture | Decent | Clean separation of concerns, good use of typed state |
| Correctness | **Poor** | Frameshift bug in 6xHis tag, validation on wrong DNA, CUA vs CTA bug |
| Biological fidelity | **Poor** | Fake backbone, meaningless RNA fold, wrong protease sites |
| Test coverage | Weak | One compound, no LLM paths, no remediation loop |
| Production readiness | Not ready | Stubs for 2 key systems, fragile LLM parsing, no retry logic |
| Spec adherence | Mixed | Pipeline shape matches, architecture diverges from spec |

The most dangerous issue is **bug #1** -- the 6xHis frameshift. Every single cassette this pipeline produces encodes a frameshifted, nonfunctional protein. If someone trusted this output and sent it to synthesis, they'd waste time and money. Bug #16 (CUA vs CTA) means a validation check is silently non-functional. Together these indicate the code was never run end-to-end against a real compound and compared to a known-good reference.
