

|  AIXBio: A LangGraph Agent Tool Suite for Biological Engineering Workflows[^1]  |
| ----- |
|   Andre van Dam **With** Apart Research  **Abstract** Biological engineering workflows — from protein sequence retrieval to expression-ready plasmid design — involve multi-step processes that demand specialized knowledge and are difficult to automate safely. We propose that LangGraph-based agent tool suites offer a promising architecture for assisting researchers with such workflows, and present AIXBio as a proof of concept. AIXBio is a two-level state graph that takes a UniProt protein accession and walks a researcher through codon optimization, expression cassette assembly, plasmid construction, and multi-point sequence validation — all orchestrated as composable graph nodes with typed state, human-in-the-loop checkpoints, and a full audit trail. The core pipeline is deterministic: an LLM is only invoked as a bounded escalation agent when rule-based remediation fails, and human review gates ensure the researcher remains in control. We demonstrate the concept on insulin (P01308) and discuss the architectural properties — composability, auditability, bounded AI agency, and testability — that make LangGraph well-suited for safety-critical biological tool suites. AIXBio is not a production system; it is a proof of concept intended to explore how agent frameworks can structure and govern multi-step bio workflows while keeping the researcher in the loop.  |

## **1\. Introduction**

Biological engineering workflows are inherently multi-step: a researcher might retrieve a protein sequence, optimize codons for a target host, design an expression cassette, assemble a plasmid, and validate the construct against multiple quality criteria — all before placing a gene synthesis order. Each step requires domain expertise, and errors compound downstream. Today, researchers either perform these steps manually using a patchwork of disconnected tools (UniProt, JCat, Benchling, SnapGene) or rely on commercial black-box services that obscure their decision logic.

Meanwhile, LLM-based agent frameworks like LangGraph [1] have demonstrated the ability to orchestrate complex, multi-step workflows with typed state management, conditional routing, human-in-the-loop interrupts, and composable subgraphs. These properties map naturally onto the requirements of biological engineering: workflows are sequential with branching logic, intermediate results must be inspectable, and human oversight is non-negotiable for safety-critical decisions.

We hypothesize that LangGraph-based agent tool suites could serve as a useful paradigm for assisting researchers with biological workflows — providing structured automation where appropriate, preserving human control at critical junctures, and maintaining a complete audit trail. To explore this hypothesis, we built AIXBio: a proof-of-concept pipeline that converts a UniProt protein accession into a validated plasmid design for *E. coli* expression.

AIXBio is not intended as a production-grade tool. Rather, it is a concrete instantiation of an architectural pattern — **deterministic-first, LLM-as-escalation, human-in-the-loop** — that we believe generalizes to a broad class of biological workflows.

**Our main contributions are:**

1. A proof-of-concept LangGraph agent tool suite that demonstrates how graph-based agent architectures can structure multi-step biological workflows with typed state, composable nodes, and human checkpoints.
2. An architectural pattern — deterministic core with bounded LLM escalation — that confines AI decision-making to a narrow, validated role while keeping the researcher in control.
3. A discussion of the properties that make agent tool suites well-suited for safety-critical bio workflows (composability, auditability, testability, bounded agency) and the open questions that remain.

## **2\. Related Work**

**Biological tools.** Codon optimization tools such as JCat [2] and the Kazusa Codon Usage Database [3] address individual steps in the expression engineering pipeline. Commercial platforms (IDT, GenScript, Thermo Fisher GeneArt) provide end-to-end services but operate as black boxes. Interactive design environments like Benchling and SnapGene offer rich interfaces but require expert operation at each step and do not compose into automated pipelines. None of these tools provide structured audit trails or programmatic human-in-the-loop checkpoints.

**AI for biology.** Recent work on AI-assisted protein engineering — ProtGPT2 [4], ESM-2 [5], and AlphaFold [6] — focuses on prediction and *de novo* design tasks upstream of expression engineering. These tools generate or evaluate protein sequences but do not address the downstream workflow of converting a chosen sequence into a clonable construct. Our work is complementary: AIXBio consumes the output of such tools (a protein sequence) and handles the subsequent engineering steps.

**Agent frameworks.** LangGraph [1] and similar frameworks (CrewAI, AutoGen) have been applied to software engineering, data analysis, and research tasks, but their application to biological engineering workflows has been limited. Bioinformatics pipelines (Nextflow, Snakemake) provide workflow orchestration but lack native support for LLM integration, human-in-the-loop interrupts, or typed state management with conditional routing. Our proof of concept explores whether LangGraph's specific combination of features — state graphs, conditional edges, interrupt/resume, and LLM tool calling — offers advantages for bio workflows.

## **3\. Methods**

### 3.1 Design Principles

We designed AIXBio around four principles that we believe are necessary for agent-assisted biological workflows:

1. **Deterministic by default.** Every step that can be performed algorithmically should be. LLMs introduce non-determinism, cost, and latency; they should only be invoked when rule-based approaches are insufficient.
2. **Bounded AI agency.** When an LLM is invoked, its scope must be precisely defined: structured inputs, typed outputs (discriminated unions), validated responses, and a maximum invocation count.
3. **Human-in-the-loop.** The researcher must be able to inspect intermediate results and approve or reject pipeline decisions at critical junctures.
4. **Full auditability.** Every decision, remediation action, and warning must be logged in the pipeline state for post-hoc review.

### 3.2 Architecture

AIXBio is implemented as a two-level LangGraph StateGraph:

**Main graph:** sequence retrieval -> human checkpoint -> fan-out to per-chain processing -> merge results -> human checkpoint -> optional structural validation.

**Chain subgraph** (runs once per protein chain): codon optimization -> cassette assembly -> plasmid assembly -> seven-point validation. On failure: deterministic remediation loop (up to 3 attempts) -> optional LLM escalation (fires at most once).

State is managed via TypedDict definitions with Annotated reducers for list accumulation fields. Nodes return partial state dicts (only the keys they modify), which enables independent testing and prevents unintended state mutations.

### 3.3 Pipeline Steps (Proof-of-Concept Implementation)

**Sequence retrieval.** Fetches the UniProt entry via REST API, extracts mature chains using annotated features (Chain, Signal, Transit peptide), and presents the result for human review. Fully deterministic.

**Codon optimization.** Greedy best-codon selection using *E. coli* K-12 usage frequencies from the Kazusa database [3], followed by a five-pass restriction site avoidance step that substitutes synonymous codons to eliminate recognition sites of the host's native enzymes (looked up via Biopython's REBASE interface). Scored by CAI [7] and GC content.

**Cassette assembly.** Concatenates: ATG start codon + 6xHis affinity tag + protease cleavage site (Enterokinase or TEV) + optimized gene + double TAA stop. Scans for N-linked glycosylation motifs and warns the researcher, since *E. coli* cannot glycosylate.

**Plasmid assembly.** Inserts the cassette into a pET-28a(+) backbone with BamHI/XhoI cloning site flanks. Generates a GenBank record with annotated features (T7 promoter, lac operator, RBS, KanR, origins of replication).

**Seven-point validation.** Each construct is evaluated against:

| Check | Threshold | Rationale |
|-------|-----------|-----------|
| GC content | 0.50–0.60 | Stability, transcription efficiency |
| CAI score | > 0.80 | Translation efficiency [7] |
| Restriction sites | 0 internal hits | Cloning compatibility |
| 5' RNA ΔG | > −10 kcal/mol | Ribosome accessibility [8] |
| Back-translation | 100% identity | Translation fidelity |
| Rare codons | 0 | Translational pausing [9] |
| Direct repeats | 0 (≥ 20 bp) | Genomic instability [10] |

*Table 1: Validation checks implemented in the proof of concept.*

**Remediation loop.** On failure, plans and applies targeted synonymous codon swaps (restriction site elimination, rare codon replacement, CAI improvement, GC adjustment). Reassembles and revalidates up to 3 times.

**LLM escalation (optional).** If deterministic remediation is exhausted, an LLM agent fires once per chain. It receives the full remediation history as structured JSON and must return one of four typed outcomes: `apply_plan` (propose swaps the greedy loop missed), `incompatible` (host cannot express this protein), `change_strategy` (change tag/vector/protease and re-run), or `give_up` (unrecoverable). Responses are validated — synonymy is checked before any swap is applied. The `escalation_used` flag prevents re-entry.

### 3.4 Repository Structure

A central goal of the proof of concept is to demonstrate that a LangGraph bio tool suite can be organized so that adding new tools and pipeline steps is straightforward — even for contributors who are not familiar with the full codebase. The repository is split into six layers, each with a single responsibility and no upward dependencies:

```
aixbio/
├── models/          Frozen dataclasses — the data contracts between steps
│   ├── protein.py       Chain, ProteinRecord
│   ├── dna.py           DNAChain, CassetteChain, CassetteElement
│   ├── plasmid.py       PlasmidChain, PlasmidRecord
│   ├── validation.py    CheckResult, ChainValidation, ValidationReport
│   ├── remediation.py   RemediationAction, PlannedFix, RemediationPlan
│   ├── escalation.py    EscalationDecision (discriminated union on `kind`)
│   ├── structure.py     StructureResult, StructureReport
│   └── audit.py         AgentDecision
├── state/           TypedDict state definitions with LangGraph reducers
│   ├── pipeline_state.py    PipelineState (main graph)
│   └── chain_state.py       ChainSubgraphState (per-chain subgraph)
├── tools/           Pure functions — no state, no side effects beyond I/O
│   ├── codon_tables.py      Codon translation, synonymous alternatives, usage table
│   ├── cai.py               Codon Adaptation Index (Biopython)
│   ├── gc.py                GC content
│   ├── restriction_sites.py REBASE enzyme lookup, recognition site scanning
│   ├── rna_fold.py          5' MFE estimation (ViennaRNA)
│   ├── repeats.py           Direct repeat detection (k-mer indexing)
│   ├── genbank.py           Plasmid GenBank file generation
│   ├── uniprot.py           UniProt REST API client
│   ├── alphafold.py         AlphaFold DB client + RMSD alignment
│   └── esmfold.py           ESMFold public API client
├── nodes/           LangGraph node functions — each reads state, returns a partial dict
│   ├── codon_optimization.py
│   ├── cassette_assembly.py
│   ├── plasmid_assembly.py
│   ├── sequence_validation.py
│   ├── remediation_agent.py     (includes apply_fixes)
│   ├── escalation_agent.py
│   ├── sequence_retrieval.py
│   ├── structural_validation.py
│   ├── human_checkpoints.py
│   ├── merge_results.py         (includes package_result, halt_pipeline, etc.)
│   └── routers.py               Conditional routing logic
├── graph/           Graph wiring only — no business logic
│   ├── main_graph.py
│   └── chain_subgraph.py
├── prompts/         LLM prompt templates (decoupled from node logic)
└── config.py        Pipeline defaults and OpenRouter LLM wrapper
```

The dependency flow is strictly downward: `graph/` imports `nodes/`, `nodes/` import `tools/` and `models/`, `tools/` import only `models/` and external libraries. This means a contributor can work on a tool without understanding the graph wiring, or add a node without modifying the models layer.

**Models** are frozen dataclasses with no logic — they define the typed contracts that flow between steps. For example, `CheckResult(name, passed, value, threshold)` is the atomic unit of validation, and `ChainValidation(id, passed, checks)` collects them per chain. Escalation outcomes use a discriminated union on a `kind` field (`"apply_plan"`, `"incompatible"`, `"change_strategy"`, `"give_up"`), so routers branch on the type tag rather than parsing free text.

**State** definitions use Python's `TypedDict` with LangGraph's `Annotated` type for list fields that accumulate across steps. An `append_log` reducer means each node returns only its *delta* (e.g., this round's remediation actions), and the framework concatenates them automatically. This prevents a common bug where a node accidentally overwrites the full history.

**Tools** are pure functions with no dependency on pipeline state. `compute_gc(dna) -> float` takes a string and returns a number; `find_restriction_sites(dna, enzymes) -> list` takes a sequence and enzyme names and returns hit positions. This makes them independently testable and reusable outside the pipeline.

**Nodes** are the glue: each is a function with the signature `f(state: ChainSubgraphState) -> dict` that reads what it needs from state, calls tools, and returns a partial state dict with only the keys it updates.

**Graph** files contain only wiring — `add_node`, `add_edge`, `add_conditional_edges` — with no business logic. The chain subgraph (`chain_subgraph.py`) registers nodes, connects them with edges, and defines conditional routing via router functions.

### 3.5 Extensibility: Adding New Tools and Nodes

The layered architecture is designed so that extending the pipeline follows a repeatable pattern. We describe the steps for two common extension scenarios.

**Adding a new validation check** (e.g., a codon pair bias score):

1. **Add a tool** in `tools/`: write a pure function, e.g., `compute_codon_pair_bias(dna: str) -> float`, in a new file `tools/codon_pair_bias.py`. The function takes a DNA string and returns a score. No pipeline imports needed — just the biology.
2. **Call it from the validation node** in `nodes/sequence_validation.py`: import the function, compute the score, and append a new `CheckResult(name="codon_pair_bias", passed=score > threshold, value=score, threshold=">0.0")` to the checks list.
3. **Optionally add a remediation strategy** in `nodes/remediation_agent.py`: add a branch that plans codon swaps when the new check fails. The remediation agent already iterates over `failed_checks` by name, so a new `elif check.name == "codon_pair_bias":` block is sufficient.

No changes to models, state, or graph wiring are required — the existing `CheckResult` dataclass and validation routing logic handle the new check automatically.

**Adding a new pipeline step** (e.g., a promoter strength prediction node):

1. **Define a model** in `models/`: create a frozen dataclass for the output, e.g., `PromoterScore(predicted_strength: float, method: str)`.
2. **Add a state field** in `state/chain_state.py`: add `promoter_score: PromoterScore | None` to `ChainSubgraphState`.
3. **Write a tool** in `tools/`: a pure function that computes the score, e.g., `predict_promoter_strength(construct_dna: str) -> float`.
4. **Write a node** in `nodes/`: a function `promoter_prediction(state: ChainSubgraphState) -> dict` that calls the tool and returns `{"promoter_score": PromoterScore(...)}`.
5. **Wire it into the graph** in `graph/chain_subgraph.py`: call `g.add_node("promoter_prediction", promoter_prediction)` and insert an edge at the desired position (e.g., between `plasmid_assembly` and `sequence_validation`).

Each layer changes independently, and the existing test patterns (`test_tools.py` for pure functions, `test_deterministic_nodes.py` for nodes in isolation) extend naturally to the new components.

This separation is what makes the tool suite a useful foundation for a broader set of bio workflows: the `tools/` layer is a growing library of reusable biological computations, the `nodes/` layer adapts them to specific pipeline steps, and the `graph/` layer composes steps into workflows that can be rearranged without rewriting the underlying logic.

### 3.6 Why LangGraph

Several LangGraph features proved particularly useful for this proof of concept:

- **StateGraph with typed state** allowed us to define explicit contracts between nodes and catch integration errors at development time rather than runtime.
- **Conditional edges and routers** naturally expressed the validation -> remediation -> escalation branching logic.
- **Send-based fan-out** handled multi-chain proteins (like insulin's A and B chains) without custom parallelization logic.
- **Interrupt/resume** provided human-in-the-loop checkpoints without additional infrastructure.
- **Annotated reducers** for list fields meant nodes could return deltas rather than full state, making the pipeline composable and testable.

These features are not unique to biology, but the combination maps well onto the structure of bio workflows where steps are sequential, branching is conditional on validation results, and human oversight is required at defined points.

## **4\. Results**

We demonstrate AIXBio on human insulin (UniProt P01308) as a representative multi-chain protein. The results below illustrate that the proof of concept is functional, not that it is comprehensive.

### 4.1 End-to-End Execution

The pipeline correctly extracts both mature insulin chains (B-chain: 30 aa, A-chain: 21 aa) from the UniProt entry, optimizes codons, assembles cassettes and plasmids, and passes all seven validation checks on the first attempt for both chains. Total execution time is under 10 seconds (excluding optional structural validation).

### 4.2 Validation Outcomes

| Check | B-chain | A-chain | Status |
|-------|---------|---------|--------|
| GC content | 0.56 | 0.54 | Pass |
| CAI score | 0.87 | 0.85 | Pass |
| Restriction sites | 0 | 0 | Pass |
| 5' RNA ΔG | −3.2 kcal/mol | −2.8 kcal/mol | Pass |
| Back-translation | 100% | 100% | Pass |
| Rare codons | 0 | 0 | Pass |
| Direct repeats | 0 | 0 | Pass |

*Table 2: Validation results for insulin. Both chains pass all checks without remediation.*

### 4.3 Architectural Observations

Building the proof of concept surfaced several properties of the LangGraph-based approach that we believe generalize:

**Composability and extensibility.** The six-layer architecture (models -> state -> tools -> nodes -> graph -> prompts) means each component can be developed, tested, and swapped independently. Adding a new validation check requires only a pure function in `tools/` and a few lines in the existing validation node — no changes to models, state, or graph wiring. Adding an entirely new pipeline step follows a repeatable five-step pattern (model, state field, tool, node, graph edge) where each layer changes independently. The test suite includes deterministic node tests that run without any LLM or network calls, and the same patterns extend naturally to new components.

**Auditability.** The pipeline state at any point contains a complete record of decisions made, remediation actions applied, and warnings emitted. For a governance use case, this means every construct can be traced back through its full decision history.

**Bounded AI agency.** The LLM escalation agent demonstrates a pattern where AI involvement is precisely scoped: it fires at most once per chain, must return one of four typed outcomes, and its proposals are validated (synonymy-checked) before execution. This is a deliberate architectural choice — the LLM assists the researcher; it does not act autonomously.

**Human-in-the-loop.** The two interrupt points (after chain extraction and after plasmid assembly) allow researchers to review and approve results before the pipeline proceeds. These are native LangGraph interrupts, not bolted-on confirmation dialogs.

### 4.4 What This Does Not Show

This proof of concept does not demonstrate: robustness across diverse proteins, benchmarking against commercial tools, wet-lab validation of produced constructs, or the escalation agent's effectiveness on real-world remediation failures. These would be necessary for a production system but are beyond the scope of exploring the architectural pattern.

## **5\. Discussion and Limitations**

### The Case for Agent Tool Suites in Biology

We believe the proof of concept illustrates a broader opportunity: LangGraph-based agent tool suites can provide structure, transparency, and governance to biological workflows that are currently either manual or opaque. The key insight is not that AIXBio automates plasmid design — commercial tools already do that — but that the agent framework provides architectural properties (typed state, conditional routing, human interrupts, audit trails) that are independently valuable for safety-critical bio workflows.

Consider the design space beyond our proof of concept. A LangGraph tool suite could assist with: CRISPR guide RNA design and off-target scoring, metabolic pathway engineering, biosafety screening of synthetic sequences, or experimental protocol planning. In each case, the pattern is similar: deterministic computation where possible, LLM assistance where judgment is needed, human approval at critical gates, and full auditability.

For AI safety specifically, this architecture creates natural governance points. An institutional biosafety committee could require that all construct designs pass through an auditable pipeline with mandatory human checkpoints — something that is difficult to enforce with ad-hoc manual workflows but straightforward with a graph-based tool suite.

### **Limitations**

- **Proof of concept only.** AIXBio has been tested on a handful of proteins. It is not validated for production use, and its biological outputs have not been wet-lab tested.
- **E. coli only.** The implementation targets a single host organism. Generalization would require new codon tables, vectors, and validation thresholds per host.
- **Simplified biology.** The pipeline does not handle many real-world challenges: codon context effects, mRNA stability beyond 5' structure, protein solubility, inclusion body risk, or disulfide bond engineering.
- **LLM escalation unvalidated.** The escalation agent's decision quality has been tested only with mocked responses. Its real-world effectiveness at resolving remediation failures is unknown.
- **Single workflow.** We demonstrate only the accession-to-plasmid workflow. The claim that the architectural pattern generalizes to other bio workflows remains a hypothesis.

### **Future Work**

- **Additional workflows.** Implement LangGraph tool suites for other bio workflows (CRISPR design, pathway engineering, biosafety screening) to test whether the architectural pattern generalizes.
- **User studies.** Evaluate whether researchers find agent-assisted workflows more efficient and less error-prone than manual approaches.
- **Governance integration.** Explore how agent tool suites can interface with institutional review processes — e.g., automatically submitting audit trails to biosafety committees.
- **Multi-host support.** Extend the proof of concept to yeast and mammalian expression systems to test composability of the node-based architecture.
- **Escalation agent evaluation.** Systematically evaluate LLM escalation quality across models, prompting strategies, and failure modes.

## **6\. Conclusion**

We presented AIXBio, a proof-of-concept LangGraph agent tool suite that walks a researcher through the multi-step process of converting a UniProt protein accession into a validated plasmid design. The contribution is not the pipeline itself — which is simplified and narrowly scoped — but the demonstration that graph-based agent architectures provide a natural fit for biological workflows that require deterministic computation, bounded AI assistance, human oversight, and full auditability.

We believe this architectural pattern — composable graph nodes, typed state, deterministic-first processing, bounded LLM escalation, and human-in-the-loop checkpoints — deserves further exploration as a paradigm for building safe, governed AI tool suites for biological research. The open questions are whether the pattern generalizes across workflows, whether researchers adopt it in practice, and how it can be integrated into institutional governance structures.

## **Code and Data**

- **Code repository**: [github.com/rht/apart\_aixbio](https://github.com/rht/apart_aixbio)
- **Data/Datasets**: Protein sequences retrieved from UniProt at runtime; no static datasets required.
- **Other artifacts**: GenBank plasmid files and validation summaries are generated in the `output/` directory for each pipeline run.

## **Author Contributions** 

A. van Dam conceived the project, designed the architecture, implemented the pipeline, and wrote this report.

## **References**

1. LangGraph: Build stateful, multi-actor applications with LLMs. LangChain, Inc. https://github.com/langchain-ai/langgraph
2. Grote, A., Hiller, K., Scheer, M., Munch, R., Nortemann, B., Hempel, D. C., & Jahn, D. (2005). JCat: a novel tool to adapt codon usage of a target gene to its potential expression host. *Nucleic Acids Research*, 33(Web Server issue), W526–W531. https://doi.org/10.1093/nar/gki376
3. Nakamura, Y., Gojobori, T., & Ikemura, T. (2000). Codon usage tabulated from international DNA sequence databases: status for the year 2000. *Nucleic Acids Research*, 28(1), 292. https://doi.org/10.1093/nar/28.1.292
4. Ferruz, N., Schmidt, S., & Hocker, B. (2022). ProtGPT2 is a deep unsupervised language model for protein design. *Nature Communications*, 13, 4348. https://doi.org/10.1038/s41467-022-32007-7
5. Lin, Z., Akin, H., Rao, R., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*, 379(6637), 1123–1130. https://doi.org/10.1126/science.ade2574
6. Jumper, J., Evans, R., Pritzel, A., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596, 583–589. https://doi.org/10.1038/s41586-021-03819-2
7. Sharp, P. M., & Li, W. H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. https://doi.org/10.1093/nar/15.3.1281
8. Kudla, G., Murray, A. W., Tollervey, D., & Plotkin, J. B. (2009). Coding-sequence determinants of gene expression in *Escherichia coli*. *Science*, 324(5924), 255–258. https://doi.org/10.1126/science.1170160
9. Kane, J. F. (1995). Effects of rare codon clusters on high-level expression of heterologous proteins in *Escherichia coli*. *Current Opinion in Biotechnology*, 6(5), 494–500. https://doi.org/10.1016/0958-1669(95)80082-4
10. Bzymek, M., & Lovett, S. T. (2001). Instability of repetitive DNA sequences: The role of replication in multiple mechanisms. *Proceedings of the National Academy of Sciences*, 98(15), 8319–8325. https://doi.org/10.1073/pnas.111008398

## **Appendix**

### A. Pipeline Architecture Diagram

```
Main Graph
├── sequence_retrieval           UniProt REST API, mature chain extraction
├── human_checkpoint             Researcher reviews extracted chains
├── fan_out_to_chains            LangGraph Send() per chain
│   └── Chain Subgraph × N
│       ├── codon_optimization   Greedy best-codon + restriction avoidance
│       ├── cassette_assembly    ATG + tag + protease + gene + stop
│       ├── plasmid_assembly     Insert into pET-28a(+), GenBank export
│       ├── sequence_validation  7 independent checks
│       └── [on failure]
│           ├── remediation_loop Synonymous swaps, up to 3 attempts
│           └── escalation_agent LLM fires once, 4 typed outcomes
├── merge_results                Aggregate per-chain outcomes
├── human_checkpoint             Researcher reviews final constructs
└── structural_validation        Optional ESMFold/AlphaFold check
```

### B. LangGraph Features Used

| Feature | Use in AIXBio | Relevance to Bio Workflows |
|---------|---------------|---------------------------|
| StateGraph | Pipeline orchestration | Explicit workflow structure |
| TypedDict + Annotated | State contracts between nodes | Prevents silent data errors |
| Conditional edges | Validation -> remediation routing | Branch on quality criteria |
| Send fan-out | Multi-chain processing | Parallelism for multi-subunit proteins |
| Interrupt/resume | Human review checkpoints | Researcher stays in the loop |
| Annotated reducers | Append-only audit logs | Immutable decision history |

### C. Key Dependencies

| Package | Purpose |
|---------|---------|
| LangGraph ≥0.2 | Graph-based pipeline orchestration with interrupts |
| Biopython ≥1.87 | Sequence translation, restriction enzymes, CAI, structure parsing |
| ViennaRNA ≥2.7.2 | RNA secondary structure MFE calculation |
| httpx ≥0.27 | Async HTTP for UniProt, ESMFold, AlphaFold APIs |

## **LLM Usage Statement**

We used Claude (Anthropic) as a pair-programming tool throughout development and to assist with drafting sections of this report. The pipeline architecture, design decisions, and all reported results were independently conceived and verified. Claude Code was used for code implementation, debugging, and iterating on the prose.


[^1]:  Research conducted at the [AIxBio Hackathon](https://apartresearch.com/sprints/aixbio-hackathon-2026-04-24-to-2026-04-26), April 2026
