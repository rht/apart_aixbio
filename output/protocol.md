# Protocol: Human Proinsulin Expression and Purification in *E. coli*

## Literature Context
Human insulin was the first recombinant therapeutic protein produced in *E. coli* (PMID: 33340315). Early work expressed the A and B chains separately as β‑galactosidase fusions, then combined and re‑oxidised to form active insulin (PMID: 85300). Modern approaches often express proinsulin (B‑C‑A) as a single polypeptide, which is subsequently converted to mature insulin by enzymatic removal of the C‑peptide (PMID: 25270715). *E. coli* remains the most cost‑effective host, but the six cysteines in proinsulin require oxidative refolding to form the correct three disulfide bonds (PMID: 38479588). A dedicated refolding step is essential for high‑yield production of bioactive insulin (PMID: 30735706).

## 1. Strains, Plasmid, and Safety Notes
- **Host strain**: *E. coli* BL21(DE3) (Novagen) – contains T7 RNA polymerase for IPTG‑inducible expression.
- **Plasmid**: pET‑28a(+) carrying the human proinsulin gene (86 aa, sequence from Uniprot P01308, residues 25–110) fused to an N‑terminal 6×His tag followed by an Enterokinase cleavage site (DDDDK). Cloning sites: *BamHI* / *XhoI*. Selection marker: kanamycin (50 µg/mL).
- **Safety**: Biosafety Level 1 (BL1). Standard aseptic technique. All waste containing IPTG or imidazole should be disposed of according to institutional guidelines.

## 2. Transformation
1. Thaw 50 µL of chemically competent BL21(DE3) cells on ice.
2. Add 1–2 µL (≈50 ng) of plasmid DNA, mix gently, incubate on ice for 30 min.
3. Heat‑shock at 42°C for 45 s, then place on ice for 2 min.
4. Add 950 µL of LB medium (no antibiotic) and incubate at 37°C, 200 rpm for 1 h.
5. Spread 100 µL onto LB agar plates containing 50 µg/mL kanamycin. Incubate at 37°C overnight.

## 3. Pre‑culture
1. Inoculate a single colony into 5 mL of LB medium + 50 µg/mL kanamycin.
2. Incubate at 37°C, 200 rpm for 12–16 h (overnight).
3. Measure OD₆₀₀; the culture should reach OD₆₀₀ ≈ 2–4.

## 4. Expression Culture and Induction
1. Inoculate 1 L of LB medium (with 50 µg/mL kanamycin) with the pre‑culture to a starting OD₆₀₀ of 0.05–0.1.
2. Grow at 37°C, 200 rpm until OD₆₀₀ reaches 0.6–0.8.
3. Induce with **0.5 mM IPTG** (final concentration).  
   *Rationale*: Predicted solubility score >0.45 (0.594) suggests high soluble expression potential, but the six cysteines (disulfide risk) favour inclusion body formation. Induction at 37°C for 4 h yields high cell mass; the protein will accumulate in inclusion bodies, which is acceptable because a refolding step is required.
4. Continue shaking at 37°C, 200 rpm for 4 h.
5. Harvest cells by centrifugation at 6,000 × g, 4°C, 15 min. Discard supernatant. Cell pellet can be stored at –80°C.

## 5. Cell Harvest and Lysis
**Denaturing lysis is required** because proinsulin forms inclusion bodies under standard induction conditions (PMID: 30735706).

1. Resuspend the cell pellet in **50 mL of denaturing lysis buffer** per 1 L culture:  
   - 50 mM Tris‑HCl, pH 8.0  
   - 300 mM NaCl  
   - 8 M urea  
   - 10 mM imidazole  
   - 1 mM β‑mercaptoethanol (optional, to reduce disulfides)
2. Stir at room temperature for 30–60 min until the suspension is homogeneous.
3. Centrifuge at 20,000 × g, 4°C, 30 min to remove insoluble debris.
4. Filter the supernatant through a 0.45 µm filter before loading onto the IMAC column.

## 6. Affinity Purification (IMAC)
All steps performed at room temperature under denaturing conditions.

1. Equilibrate a Ni‑NTA agarose column (5 mL resin per 1 L culture) with **Buffer A**: 50 mM Tris‑HCl, pH 8.0, 300 mM NaCl, 8 M urea, 10 mM imidazole.
2. Load the clarified lysate at a flow rate of 1 mL/min.
3. Wash with 10 column volumes (CV) of Buffer A.
4. Wash with 5 CV of **Buffer B**: 50 mM Tris‑HCl, pH 8.0, 300 mM NaCl, 8 M urea, 20 mM imidazole.
5. Elute with **Buffer C**: 50 mM Tris‑HCl, pH 8.0, 300 mM NaCl, 8 M urea, 250 mM imidazole. Collect 1 mL fractions.
6. Pool fractions containing proinsulin (monitor by A₂₈₀ or SDS‑PAGE). Typical yield: 10–30 mg per L culture.

## 7. Tag Removal
1. Dialyse the pooled eluate against **refolding buffer** (see Section 8) to remove imidazole and urea gradually.  
   *Alternatively*, perform tag removal after refolding to avoid premature cleavage in denaturant.
2. After refolding (Section 8), add **Enterokinase** (e.g., recombinant enterokinase, 1 U per 50 µg of fusion protein) and incubate at 25°C for 16 h.
3. Remove the cleaved His‑tag and uncleaved protein by passing through a second Ni‑NTA column (equilibrated with native buffer). The flow‑through contains tag‑free proinsulin.
4. *Note*: To obtain mature insulin, the C‑peptide must be removed by subsequent treatment with trypsin and carboxypeptidase B (not covered in this protocol).

## 8. Inclusion Body Refolding
Because **disulfide_risk = True** (6 cysteines), oxidative refolding is mandatory to form the three native disulfide bonds (B7‑A7, B19‑A20, A6‑A11).

### a. Solubilisation of inclusion bodies (if not already in denaturing buffer)
- Resuspend the IMAC‑purified proinsulin in **solubilisation buffer**: 50 mM Tris‑HCl, pH 8.0, 6 M guanidine‑HCl, 10 mM DTT. Incubate at 37°C for 1 h.
- Centrifuge at 20,000 × g, 30 min to remove any insoluble material.

### b. Refolding by dilution
1. Prepare **refolding buffer** (freshly made, degassed):  
   - 50 mM Tris‑HCl, pH 8.5  
   - 0.5 M L‑arginine (suppresses aggregation)  
   - 5 mM reduced glutathione (GSH)  
   - 0.5 mM oxidized glutathione (GSSG)  
   - 1 mM EDTA  
   - 0.1 mM PMSF (protease inhibitor)
2. Slowly add the denatured protein dropwise to the refolding buffer at 4°C with gentle stirring, to a final protein concentration of **0.1–0.2 mg/mL**.
3. Incubate at 4°C for 24–48 h without stirring.
4. Dialyse against **50 mM Tris‑HCl, pH 8.0, 150 mM NaCl** to remove arginine and glutathione.
5. Centrifuge at 20,000 × g, 30 min to remove any precipitate. The supernatant contains refolded proinsulin.

### c. Quality check of refolding
- Measure thiol content using Ellman’s assay; free thiols should be <0.1 per molecule.
- Analyse by non‑reducing SDS‑PAGE: correctly folded proinsulin migrates slightly faster than the reduced form.

## 9. Quality Control
- **SDS‑PAGE** (reducing and non‑reducing): expected molecular weight of proinsulin ≈ 9.4 kDa (without tag) or ≈ 11.5 kDa (with His‑tag).
- **Western blot**: anti‑insulin antibody (recognises B‑chain) or anti‑His antibody.
- **Mass spectrometry**: confirm intact mass (expected 9,398 Da for tag‑free proinsulin).
- **Insulin activity assay**: after conversion to mature insulin, measure glucose uptake in adipocytes or use ELISA (PMID: 30735706).

## 10. Expected Yields and Troubleshooting Tips
- **Expected yield**: 10–30 mg of purified proinsulin per litre of culture (after IMAC and refolding). Final recovery after tag removal and conversion may be 5–15 mg/L.
- **Low soluble expression**: If soluble proinsulin is desired (e.g., for periplasmic expression), consider using a different strain (e.g., *E. coli* 20, PMID: 30735706) or lower induction temperature (25°C, 0.1 mM IPTG). However, refolding will still be needed for disulfide formation.
- **Aggregation during refolding**: Reduce protein concentration to ≤0.1 mg/mL; increase L‑arginine to 1 M; adjust pH to 8.5–9.0; or use a redox shuffling system (e.g., 2 mM GSH/0.4 mM GSSG).
- **Incomplete disulfide formation**: Increase incubation time to 72 h; add a catalytic amount of protein disulfide isomerase (PDI) if available.
- **Low IMAC binding**: Ensure urea concentration is ≥6 M; imidazole in lysis buffer should be ≤10 mM; regenerate Ni‑NTA resin after each use.

**Note**: For therapeutic use, additional polishing steps (ion‑exchange, SEC) and endotoxin removal are required. Conversion to mature insulin is performed by sequential digestion with trypsin and carbox