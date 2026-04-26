


|  BioClaw: An Agentic Tool Suite for Biological Engineering Workflows[^1]  |
| ----- |
|   Drew **With** Apart Research
**Abstract** We present BioClaw, a modular agent framework for biological engineering workflows. We demonstrate it by converting human insulin (UniProt P01308) into validated, expression-ready DNA constructs for *E. coli*. Starting from a protein accession, the pipeline retrieves the sequence, optimizes codons, assembles expression cassettes, builds annotated constructs, and runs seven-point validation. The insulin construct passes all checks after one remediation round. BioClaw is built on LangGraph as a two-level state graph with composable nodes, typed state, human-in-the-loop checkpoints, and a full audit trail. The core pipeline is deterministic: an LLM is only invoked as a bounded escalation agent when rule-based remediation fails. New validation checks, pipeline steps, or expression hosts can be added through a repeatable extension pattern without restructuring the graph. We discuss the properties that make this approach well-suited for safety-critical biological workflows and show that composability, auditability, bounded AI agency, and testability emerge naturally from the graph-based architecture.

## **1. Introduction**

Biological engineering workflows are inherently multi-step. A researcher might retrieve a protein sequence, optimize codons for a target host, design an expression cassette, assemble a plasmid, and validate the construct against multiple quality criteria, all before placing a gene synthesis order. Each step requires domain expertise, and errors compound downstream. Today, researchers either perform these steps manually using disconnected tools (UniProt, JCat, Benchling, SnapGene) or rely on commercial black-box services that obscure their decision logic.

We hypothesize that LangGraph-based [1] agent tool suites could serve as a useful paradigm for assisting researchers with biological workflows: structured automation where appropriate, human control at critical junctures, and a complete audit trail. To explore this hypothesis, we built BioClaw: a proof-of-concept pipeline that converts a UniProt protein accession into a validated plasmid design for *E. coli* expression.

BioClaw instantiates one architectural pattern, **deterministic-first, LLM-as-escalation, human-in-the-loop**, that we believe generalizes to a broad class of biological workflows.

**Our main contributions are:**

1. A LangGraph agent tool suite that demonstrates the architecture on a concrete biological workflow (UniProt accession to validated plasmid design).
2. An architectural pattern (deterministic-first, LLM-as-escalation, human-in-the-loop) that confines AI decision-making to a narrow, validated role while keeping the researcher in control.
3. A discussion of properties and open questions for agent tool suites in safety-critical biological workflows.

## **2. Related Work**

**Biological tools.** Codon optimization tools such as JCat [2] and the Kazusa Codon Usage Database [3] address individual steps in the expression engineering pipeline. Commercial platforms (IDT, GenScript, Thermo Fisher GeneArt) provide end-to-end services but operate as black boxes. Interactive design environments like Benchling and SnapGene offer rich interfaces but require expert operation at each step and do not compose into automated pipelines. None of these tools provide structured audit trails or programmatic human-in-the-loop checkpoints.

**AI for biology.** Recent work on AI-assisted protein engineering (ProtGPT2 [4], ESM-2 [5]) focuses on prediction and *de novo* design tasks upstream of expression engineering. These tools generate or evaluate protein sequences but do not address the downstream workflow of converting a chosen sequence into a clonable construct. Our work is complementary: BioClaw consumes the output of such tools (a protein sequence) and handles the subsequent engineering steps. For structural validation of the final DNA constructs, we use Evo2 [6], a genomic foundation model that provides per-nucleotide log-likelihood scores as a proxy for sequence plausibility.

## **3. Methods**

### 3.1 Design Principles

We designed BioClaw around four principles:

1. **Deterministic by default.** Every step that can be performed algorithmically should be. LLMs introduce non-determinism, cost, and latency. They should only be invoked when rule-based approaches are insufficient.
2. **Bounded AI agency.** When an LLM is invoked, its scope must be precisely defined: structured inputs, typed outputs (discriminated unions), validated responses, and a maximum invocation count.
3. **Human-in-the-loop.** The researcher must be able to inspect intermediate results and approve or reject pipeline decisions at critical junctures.
4. **Full auditability.** Every decision, remediation action, and warning must be logged in the pipeline state for post-hoc review.

### 3.2 Architecture

BioClaw is implemented as a two-level LangGraph StateGraph:

**Main graph:** sequence retrieval, then a human checkpoint, then fan-out to per-chain processing, then merge results, then a second human checkpoint, then optional structural validation.

**Chain subgraph** (runs once per protein chain): codon optimization, cassette assembly, plasmid assembly, and seven-point validation. On failure, a deterministic remediation loop revalidates up to 3 times. If that fails and escalation is enabled, an LLM agent fires once with typed outcomes (see Section 3.3).

State is managed via TypedDict definitions with Annotated reducers for list accumulation fields. Nodes return partial state dicts (only the keys they modify). This enables independent testing and prevents unintended state mutations.

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

*Table 1: Validation checks.*

**Remediation loop.** On failure, plans and applies targeted synonymous codon swaps (restriction site elimination, rare codon replacement, CAI improvement, GC adjustment). Reassembles and revalidates up to 3 times.

**LLM escalation (optional).** If deterministic remediation is exhausted, an LLM agent fires once per chain. It receives the full remediation history as structured JSON and must return one of four typed outcomes: `apply_plan` (propose swaps the greedy loop missed), `incompatible` (host cannot express this protein), `change_strategy` (change tag/vector/protease and re-run), or `give_up` (unrecoverable). Responses are validated before any swap is applied. The `escalation_used` flag prevents re-entry.

**Structural validation (optional).** When enabled, the pipeline scores each final DNA construct using Evo2-1B-base [6], a genomic foundation model accessed via the BioLM API. Evo2 returns per-nucleotide log-likelihoods; the pipeline computes a mean log-likelihood for the construct as a measure of sequence plausibility. This provides a learned, structure-aware quality signal complementary to the rule-based seven-point validation.

## **4. Results**

We demonstrate BioClaw on human insulin (UniProt P01308).

### 4.1 End-to-End Execution

The pipeline correctly extracts the mature insulin sequence (86 aa, after stripping the 24-residue signal peptide) from the UniProt entry, optimizes codons, assembles the cassette and plasmid, and passes all seven validation checks after one remediation round. Total execution time is under 10 seconds (excluding optional structural validation).

### 4.2 Validation Outcomes

| Check | Insulin_P01308 | Status |
|-------|----------------|--------|
| GC content | 0.54 | Pass |
| CAI score | 0.89 | Pass |
| Restriction sites | 0 | Pass |
| 5' RNA ΔG | -0.7 kcal/mol | Pass |
| Back-translation | 100% | Pass |
| Rare codons | 0 | Pass |
| Direct repeats | 0 | Pass |

*Table 2: Validation results for insulin (single mature chain, 86 aa). All checks pass after one remediation round.*

### 4.3 What This Does Not Show

These results cover a single protein. See Section 5 for a full discussion of limitations.

## **5. Discussion and Limitations**

### The Case for Agent Tool Suites in Biology

The situation in biological engineering resembles software development before CI/CD: repeatable, rule-based work is executed by hand and errors compound silently downstream. This matters most for biosafety screening, where an institutional biosafety committee could require that every construct pass through an auditable pipeline with mandatory human checkpoints, something difficult to enforce with ad-hoc manual workflows. The same pattern extends to CRISPR guide design, pathway engineering, or any workflow where reproducibility and oversight are non-negotiable.

### **Limitations**

- **Proof of concept only.** BioClaw has been tested on a handful of proteins. It is not validated for production use, and its biological outputs have not been wet-lab tested.
- **E. coli only.** The implementation targets a single host organism. Generalization would require new codon tables, vectors, and validation thresholds per host.
- **Simplified biology.** The pipeline does not handle many real-world challenges: codon context effects, mRNA stability beyond 5' structure, or disulfide bond engineering. The solubility predictor uses a hand-tuned composition heuristic (GRAVY, instability index, charged fraction, pI) with unvalidated weights, not a trained model; replacing it with a validated predictor such as NetSolP or Protein-Sol is straightforward future work.
- **LLM escalation unvalidated.** The escalation agent's decision quality has been tested only with mocked responses. Its real-world effectiveness at resolving remediation failures is unknown.
- **Single workflow.** We demonstrate only the accession-to-plasmid workflow. The claim that the pattern generalizes to other bio workflows remains a hypothesis.

### **Future Work**

- **Additional workflows.** Implement tool suites for other bio workflows (CRISPR design, pathway engineering, biosafety screening).
- **User studies.** Evaluate whether researchers find agent-assisted workflows more efficient and less error-prone than manual approaches.
- **Governance integration.** Explore how agent tool suites can interface with institutional review processes, for example by automatically submitting audit trails to biosafety committees.
- **Multi-host support.** Extend the proof of concept to yeast and mammalian expression systems to test composability of the node-based architecture.
- **Escalation agent evaluation.** Systematically evaluate LLM escalation quality across models, prompting strategies, and failure modes.

## **6. Conclusion**

We presented BioClaw, a LangGraph agent tool suite that walks a researcher through converting a UniProt protein accession into a validated plasmid design. The main contribution is the demonstration that graph-based agent architectures provide a natural fit for biological workflows. The open questions are whether the pattern generalizes across workflows, whether researchers adopt it in practice, and how it can be integrated into institutional governance structures.

## **Code and Data**

- **Code repository**: [github.com/rht/apart\_aixbio](https://github.com/rht/apart_aixbio)
- **Data/Datasets**: Protein sequences retrieved from UniProt at runtime. No static datasets required.
- **Other artifacts**: GenBank plasmid files and validation summaries are generated in the `output/` directory for each pipeline run.

## **Author Contributions**

A. van Dam conceived the project, designed the architecture, implemented the pipeline, and wrote this report.

## **References**

1. LangGraph. LangGraph. https://www.langchain.com/langgraph
2. Grote A, Hiller K, Scheer M, Munch R, Nortemann B, Hempel DC, Jahn D. JCat: a novel tool to adapt codon usage of a target gene to its potential expression host. Nucleic Acids Research [Internet]. 2005 Jun 26;33(Web Server):W526–W531. Available from: https://doi.org/10.1093/nar/gki376
3. Nakamura Y, Gojobori T, Ikemura T. Codon usage tabulated from international DNA sequence databases: status for the year 2000. Nucleic Acids Res. 2000 Jan 1;28(1):292. doi: 10.1093/nar/28.1.292. PMID: 10592250; PMCID: PMC102460.
4. Ferruz N, Schmidt S, Höcker B. ProtGPT2 is a deep unsupervised language model for protein design. Nature Communications [Internet]. 2022 Jul 27;13(1):4348. Available from: https://doi.org/10.1038/s41467-022-32007-7
5. Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, Smetanin N, Verkuil R, Kabeli O, Shmueli Y, Costa ADS, Fazel-Zarandi M, Sercu T, Candido S, Rives A. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science [Internet]. 2023 Mar 16;379(6637):1123–1130. Available from: https://doi.org/10.1126/science.ade2574
6. Arc Institute. Evo2. 2025. Available from: https://arcinstitute.org/news/evo2
7. Sharp PM, Li WH. The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research [Internet]. 1987 Jan 1;15(3):1281–1295. Available from: https://doi.org/10.1093/nar/15.3.1281
8. Kudla G, Murray AW, Tollervey D, Plotkin JB. Coding-Sequence Determinants of Gene Expression in Escherichia coli. Science [Internet]. 2009 Apr 9;324(5924):255–258. Available from: https://doi.org/10.1126/science.1170160
9. Kane JF. Effects of rare codon clusters on high-level expression of heterologous proteins in Escherichia coli. Current Opinion in Biotechnology [Internet]. 1995 Jan 1;6(5):494–500. Available from: https://doi.org/10.1016/0958-1669(95)80082-4
10. Bzymek M, Lovett ST. Instability of repetitive DNA sequences: The role of replication in multiple mechanisms. Proceedings of the National Academy of Sciences [Internet]. 2001 Jul 17;98(15):8319–8325. Available from: https://doi.org/10.1073/pnas.111008398

## **Appendix**

## **LLM Usage Statement**

We used Claude (Anthropic) as a pair-programming tool throughout development and to assist with drafting sections of this report. The pipeline architecture, design decisions, and all reported results were independently conceived and verified. Claude Code was used for code implementation, debugging, and iterating on the prose.


[^1]:  Research conducted at the [AIxBio Hackathon](https://apartresearch.com/sprints/aixbio-hackathon-2026-04-24-to-2026-04-26), April 2026
