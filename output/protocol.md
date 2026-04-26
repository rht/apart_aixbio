# Protocol: Human Proinsulin Expression and Purification in *E. coli*

## Literature Context
Human insulin is a therapeutic protein historically produced in *E. coli* via separate expression of A and B chains (PMID:85300) or as a single-chain proinsulin precursor (PMID:33340315). Modern production often uses *E. coli* inclusion bodies followed by oxidative refolding to form the three correct disulfide bonds (PMID:25270715). The solubility score (0.594) predicts a soluble fraction, but the presence of six cysteines (disulfide risk = True) necessitates a dedicated refolding step to achieve native disulfide pairing (PMID:30735706). Denaturing lysis and glutathione‑based redox refolding are recommended (PMID:38479588).

## 1. Strains, Plasmid, and Safety Notes
- **Host strain**: *E. coli* BL21(DE3) (Novagen) – carries T7 RNA polymerase under lacUV5 control.
- **Plasmid**: pET‑28a(+) encoding human proinsulin (86 aa) with an N‑terminal 6×His tag followed by an enterokinase (EK) cleavage site (DDDDK↓). Cloning sites: *BamHI* / *XhoI*.
- **Antibiotic**: Kanamycin (50 µg/mL) for plasmid maintenance.
- **Safety**: *E. coli* BL21(DE3) is BSL‑1. Guanidine hydrochloride and urea are irritants; work in a fume hood. β‑mercaptoethanol (BME) is toxic – handle in a chemical hood.

## 2. Transformation
1. Thaw 50 µL of chemically competent BL21(DE3) cells on ice.
2. Add 1 µL (10–50 ng) of pET‑28a‑proinsulin plasmid. Mix gently.
3. Incubate on ice for 30 min.
4. Heat‑shock at 42°C for 45 s. Place back on ice for 2 min.
5. Add 950 µL of LB medium (no antibiotic). Incubate at 37°C, 220 rpm for 1 h.
6. Spread 100 µL on LB‑agar plates containing 50 µg/mL kanamycin. Incubate at 37°C overnight.

## 3. Pre‑culture
1. Inoculate a single colony into 5 mL of LB + 50 µg/mL kanamycin.
2. Grow at 37°C, 220 rpm for 12–16 h (OD₆₀₀ ≈ 2–3).

## 4. Expression Culture and Induction
1. Inoculate 500 mL of LB + kanamycin (50 µg/mL) with pre‑culture to an initial OD₆₀₀ of 0.05–0.1.
2. Grow at 37°C, 220 rpm until OD₆₀₀ reaches 0.6–0.8.
3. **Induction**: Add IPTG to a final concentration of **0.5 mM**.
4. **Temperature**: 37°C for 4 h (high‑yield inclusion body formation).  
   *Alternative*: 25°C overnight if soluble expression is desired, but refolding will still be required.
5. Continue shaking at 220 rpm for the induction period.

## 5. Cell Harvest and Lysis
1. Harvest cells by centrifugation at 6,000 × g, 4°C, 15 min.
2. Discard supernatant. Cell pellet can be stored at –80°C.
3. **Denaturing lysis buffer**: 8 M urea, 50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 10 mM imidazole, 1 mM BME. Add fresh BME.
4. Resuspend pellet in 20 mL lysis buffer per gram wet weight.
5. Lyse by sonication on ice (6 × 30 s pulses, 50% amplitude, 2 min cooling between cycles).
6. Centrifuge at 20,000 × g, 4°C, 30 min. Retain supernatant (denatured lysate).

## 6. Affinity Purification (IMAC)
1. Equilibrate 5 mL Ni‑NTA agarose (Qiagen) with denaturing lysis buffer.
2. Load clarified lysate onto column at 1 mL/min (gravity or peristaltic pump).
3. Wash with 10 column volumes (CV) of wash buffer: 8 M urea, 50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 20 mM imidazole, 1 mM BME.
4. Elute with 5 CV of elution buffer: 8 M urea, 50 mM Tris‑HCl pH 8.0, 300 mM NaCl, 250 mM imidazole, 1 mM BME. Collect 1 mL fractions.
5. Pool fractions containing proinsulin (SDS‑PAGE, ~9.5 kDa for proinsulin alone; with His‑EK tag ~12 kDa).

## 7. Tag Removal (Enterokinase Cleavage)
*Perform after refolding to avoid interference with disulfide formation.*
1. Dialyze refolded proinsulin (from Section 8) into EK cleavage buffer: 50 mM Tris‑HCl pH 8.0, 50 mM NaCl, 2 mM CaCl₂.
2. Add recombinant enterokinase (e.g., Novagen EKMax) at 1 U per 50 µg protein.
3. Incubate at 25°C for 16 h.
4. Remove uncleaved protein and EK by passing through a second Ni‑NTA column (flow‑through contains tag‑free proinsulin). Wash with 2 CV of cleavage buffer.
5. Concentrate using a 3 kDa MWCO centrifugal filter.

## 8. Inclusion Body Refolding
*Required because disulfide_risk = True (6 cysteines).*

### a. Denaturation and Reduction
- Adjust pooled IMAC eluate to 10 mg/mL protein (Bradford assay).
- Add DTT to 10 mM. Incubate at 37°C for 1 h to fully reduce cysteines.

### b. Refolding by Dialysis
- Prepare refolding buffer (fresh): 50 mM Tris‑HCl pH 8.5, 0.5 M L‑arginine, 1 mM EDTA, 3 mM reduced glutathione (GSH), 0.3 mM oxidized glutathione (GSSG).  
  *Redox ratio 10:1 (GSH:GSSG) favours correct disulfide formation (PMID:25270715).*
- Dialyze the reduced protein against 100 volumes of refolding buffer at 4°C for 24 h, with one buffer change after 12 h.
- After dialysis, centrifuge at 20,000 × g, 4°C, 30 min to remove aggregates.

### c. Oxidative Refolding Check
- Monitor by reverse‑phase HPLC or non‑reducing SDS‑PAGE (shift to more compact band indicates disulfide formation).

## 9. Quality Control
1. **SDS‑PAGE** (reducing and non‑reducing): Expected bands – reduced: ~12 kDa (His‑proinsulin), ~9.5 kDa (after tag removal). Non‑reducing: more compact band if disulfides formed.
2. **Western blot**: Anti‑6×His antibody (before tag removal) or anti‑insulin antibody.
3. **Mass spectrometry**: Confirm intact mass (proinsulin: 9,398 Da; after tag removal: ~5,808 Da for A‑B chains? Actually proinsulin is 86 aa, ~9.4 kDa; after C‑peptide removal by EK? Wait – EK cleaves after DDDDK, not at C‑peptide junctions. The tag is removed, but proinsulin remains as a single chain. For active insulin, further processing with trypsin and carboxypeptidase B is needed (not covered here). The protocol yields proinsulin. If active insulin is required, refer to PMID:85300 for separate chain expression and combination.)
4. **Insulin activity**: Optional – ELISA or cell‑based glucose uptake assay.

## 10. Expected Yields and Troubleshooting Tips
- **Expected yield**: 10–50 mg of purified proinsulin per liter of culture (typical for inclusion body route).
- **Low solubility after refolding**: Increase L‑arginine to 1 M or adjust GSH:GSSG ratio (e.g., 5:1 to 20:1).
- **No disulfide formation**: Ensure complete reduction before refolding; add 0.1 mM DTT to refolding buffer if needed.
- **Enterokinase cleavage inefficient**: Extend incubation to 24 h or add fresh enzyme after 12 h.
- **Co‑expression of chaperones**: Not required because inclusion body route is used; if soluble expression is attempted, co‑express GroEL/GroES (PMID:38479588) and induce at 25°C with 0.1 mM IPTG.

**Note**: For production of active human insulin (A and B chains with correct inter‑chain disulfides), express A and B chains separately, refold each, then combine in a 1:1 molar ratio and perform oxidative refolding as described in PMID:85300.