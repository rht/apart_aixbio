# Protocol: Insulin (Proinsulin) Expression and Purification in *E. coli*

## Literature Context
Human insulin is a therapeutic protein historically produced in *E. coli* via separate expression of A and B chains followed by chemical combination (PMID:85300). Modern approaches often express proinsulin as a single chain, requiring oxidative refolding to form three disulfide bonds (Cys⁶–Cys¹¹, Cys⁷–Cys⁷, Cys²⁰–Cys¹⁹). The *E. coli* host is cost‑effective but lacks the eukaryotic machinery for disulfide bond formation; therefore, denaturing purification and controlled refolding are essential (PMID:38479588, PMID:25270715). The pET‑28a(+) vector with an N‑terminal 6×His tag and enterokinase cleavage site enables IMAC purification and tag removal. The predicted solubility score (0.594) suggests good intrinsic solubility, but the presence of six cysteines (disulfide risk = True) mandates an inclusion body refolding strategy to obtain correctly folded, bioactive insulin (PMID:30735706).

## 1. Strains, Plasmid, and Safety Notes
- **Host strain**: *E. coli* BL21(DE3) (or derivative, e.g., BL21‑CodonPlus‑RIL for rare codon supplementation if needed; here CAI = 0.89, so not required).
- **Plasmid**: pET‑28a(+) encoding human proinsulin (86 aa) with N‑terminal 6×His tag followed by an enterokinase recognition site (DDDDK). Cloning sites: *BamHI* (5′) and *XhoI* (3′). Kanamycin resistance (50 µg/mL).
- **Safety**: *E. coli* BL21(DE3) is BSL‑1. Use standard aseptic technique. Denaturing agents (urea, guanidine‑HCl) and reducing agents (DTT, β‑mercaptoethanol) are irritants; work in a fume hood.

## 2. Transformation
1. Thaw 50 µL of chemically competent BL21(DE3) cells on ice.
2. Add 1–2 µL (≈50 ng) of plasmid DNA, mix gently, incubate on ice for 30 min.
3. Heat‑shock at 42°C for 45 s, then place on ice for 2 min.
4. Add 950 µL of LB medium (no antibiotic) and incubate at 37°C, 220 rpm for 1 h.
5. Spread 100 µL on LB agar plates containing 50 µg/mL kanamycin. Incubate at 37°C overnight (16–18 h).

## 3. Pre‑culture
1. Inoculate a single colony into 5 mL of LB medium + 50 µg/mL kanamycin in a 15 mL tube.
2. Incubate at 37°C, 220 rpm for 12–16 h (overnight) until OD₆₀₀ ≈ 2–3.

## 4. Expression Culture and Induction
- **Induction conditions**: Because the solubility score is >0.45 but disulfide risk is high, we aim for high‑level expression as inclusion bodies. Use 1.0 mM IPTG at 37°C for 4 h (PMID:30735706).
1. Inoculate 1 L of LB medium (or TB for higher density) + 50 µg/mL kanamycin with 10 mL of pre‑culture (1% v/v).
2. Grow at 37°C, 220 rpm until OD₆₀₀ reaches 0.6–0.8.
3. Add IPTG to a final concentration of 1.0 mM.
4. Continue incubation at 37°C, 220 rpm for 4 h.
5. Harvest cells by centrifugation at 6,000 × g, 4°C, 15 min. Discard supernatant. Cell pellet can be stored at –80°C.

## 5. Cell Harvest and Lysis
- **Denaturing lysis** is required because the protein will be in inclusion bodies (warnings).
1. Resuspend cell pellet in 30 mL of **Lysis Buffer** (50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 8 M urea, 10 mM imidazole, 5 mM β‑mercaptoethanol). Use a homogeniser or magnetic stirrer.
2. Stir at room temperature for 1 h to fully solubilise inclusion bodies.
3. Centrifuge at 20,000 × g, 4°C, 30 min. Retain the supernatant (clarified lysate). Filter through a 0.45 µm syringe filter.

## 6. Affinity Purification (IMAC)
- All steps at room temperature under denaturing conditions.
1. Equilibrate 5 mL of Ni‑NTA agarose resin (e.g., Qiagen) with **Equilibration Buffer** (50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 8 M urea, 10 mM imidazole).
2. Load clarified lysate onto the column at a flow rate of 1 mL/min.
3. Wash with 10 column volumes (CV) of **Wash Buffer** (50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 8 M urea, 20 mM imidazole).
4. Elute with **Elution Buffer** (50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 8 M urea, 250 mM imidazole). Collect 1 mL fractions. Monitor A₂₈₀.
5. Pool fractions containing proinsulin (typically fractions 2–5). Check purity by SDS‑PAGE (Coomassie stain).

## 7. Tag Removal
- The 6×His tag is removed by enterokinase cleavage after refolding (to avoid interference). Perform after refolding (Section 8). Alternatively, cleave under denaturing conditions if the tag interferes with refolding; here we recommend post‑refolding.
1. After refolding (Section 8), dialyse the refolded protein into **Cleavage Buffer** (50 mM Tris‑HCl pH 8.0, 100 mM NaCl, 2 mM CaCl₂).
2. Add enterokinase (e.g., New England Biolabs) at 1 U per 50 µg of fusion protein. Incubate at 25°C for 16 h.
3. Remove uncleaved protein and enterokinase by passing through a second Ni‑NTA column (flow‑through contains cleaved proinsulin). Wash with 2 CV of cleavage buffer.
4. Concentrate the flow‑through using a 3 kDa MWCO centrifugal filter.

## 8. Inclusion Body Refolding
- **Rationale**: Six cysteines require oxidative refolding to form the three native disulfide bonds (CysA6–CysA11, CysA7–CysB7, CysA20–CysB19). The refolding protocol is adapted from classical insulin refolding (PMID:85300, PMID:25270715).
1. **Prepare refolding buffer**: 50 mM Tris‑HCl pH 8.5, 0.5 M L‑arginine (to suppress aggregation), 2 mM reduced glutathione (GSH), 0.2 mM oxidised glutathione (GSSG) (redox ratio 10:1), 1 mM EDTA.
2. Dilute the purified denatured proinsulin (in 8 M urea) dropwise into refolding buffer at 4°C with gentle stirring to a final protein concentration of 0.1–0.2 mg/mL. Avoid exceeding 0.2 mg/mL to minimise aggregation.
3. Incubate at 4°C for 24–48 h without stirring.
4. Monitor refolding by reverse‑phase HPLC or by measuring thiol groups (Ellman’s assay). The formation of disulfide bonds can be confirmed by a shift in retention time or by mass spectrometry.
5. After refolding, dialyse against 20 mM Tris‑HCl pH 8.0, 100 mM NaCl to remove arginine and glutathione.
6. Centrifuge at 20,000 × g, 4°C, 30 min to remove any precipitate. Retain supernatant.

## 9. Quality Control
1. **SDS‑PAGE** (reducing and non‑reducing): Compare mobility. Correctly folded proinsulin should show a slightly faster migration under non‑reducing conditions due to compact structure.
2. **Mass spectrometry**: Expected mass of proinsulin (86 aa) ≈ 9.4 kDa; after tag removal ≈ 8.8 kDa. Confirm intact mass and disulfide bond formation (mass shift of –6 Da for three disulfides).
3. **Reverse‑phase HPLC**: Use a C18 column with a gradient of 20–60% acetonitrile in 0.1% TFA. Native proinsulin elutes earlier than the reduced form.
4. **Bioactivity assay** (optional): Insulin receptor binding or glucose uptake assay in adipocytes.

## 10. Expected Yields and Troubleshooting Tips
- **Expected yield**: 10–30 mg of purified proinsulin per litre of culture (after refolding). Final recovery after tag cleavage and polishing may be 5–15 mg/L.
- **Low solubility after refolding**: Reduce protein concentration during refolding to ≤0.1 mg/mL; increase L‑arginine to 1 M; adjust redox ratio (try GSH:GSSG = 5:1 to 20:1).
- **Incomplete disulfide formation**: Extend refolding time to 72 h; add a redox shuffling step (e.g., 0.5 mM DTT for 1 h then removal).
- **Tag cleavage inefficiency**: Ensure Ca²⁺ is present; increase enterokinase amount or incubation time; remove imidazole completely by dialysis.
- **Aggregation during IMAC**: Keep urea concentration at 8 M throughout; do not reduce imidazole below 10 mM in wash steps.
- **Co‑expression of chaperones**: Not required here because we are using inclusion body route; if soluble expression is desired, co‑express GroEL/GroES (PMID:38479588), but for disulfide‑rich proteins the refolding route is more reliable.