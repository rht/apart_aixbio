# Protocol: Human Proinsulin Expression and Purification in *E. coli*

## Literature Context
Human insulin was the first recombinant therapeutic protein produced in *E. coli* (PMID:33340315). Early work expressed the A and B chains separately as β-galactosidase fusions, then combined and re‑oxidised to form active insulin (PMID:85300). Modern processes typically express proinsulin in *E. coli* as inclusion bodies, followed by denaturing purification and oxidative refolding (PMID:25270715). The target solubility score (0.594) predicts moderate-to-high solubility, but the presence of six cysteines (disulfide risk = True) mandates oxidative refolding to obtain the correct three disulfide bonds (PMID:25270715, PMID:38479588). A new *E. coli* strain engineered for insulin production has been described (PMID:30735706), but we use the standard BL21(DE3)/pET‑28a system.

## 1. Strains, Plasmid, and Safety Notes
- **Host strain**: *E. coli* BL21(DE3) (Novagen) – contains T7 RNA polymerase under lacUV5 control.
- **Plasmid**: pET‑28a(+) encoding the 86‑aa human proinsulin (UniProt P01308, residues 25–110) with an N‑terminal 6×His‑tag followed by an Enterokinase (EK) recognition site (DDDDK). Cloning sites: *BamHI* (5′) and *XhoI* (3′).
- **Resistance marker**: Kanamycin (50 µg/mL final).
- **Safety**: Standard BSL‑1 practices. Use PPE when handling IPTG and guanidine hydrochloride. Ni‑NTA resin is an irritant.

## 2. Transformation
1. Thaw 50 µL of chemically competent BL21(DE3) cells on ice.
2. Add 1–5 ng (1 µL) of plasmid DNA. Mix gently, incubate on ice 30 min.
3. Heat‑shock at 42 °C for 45 s, then place on ice 2 min.
4. Add 450 µL LB medium (no antibiotic). Incubate 1 h at 37 °C, 200 rpm.
5. Spread 50–100 µL on LB agar containing 50 µg/mL kanamycin. Incubate 16 h at 37 °C.
6. Pick a single colony for pre‑culture.

## 3. Pre‑culture
1. Inoculate 20 mL LB (50 µg/mL kanamycin) in a 100 mL flask with one colony.
2. Incubate at 37 °C, 200 rpm for 12–16 h (overnight) until OD₆₀₀ ≈ 2.0–3.0.
3. Use this pre‑culture to inoculate the expression culture at 1:100 dilution.

## 4. Expression Culture and Induction
- *Solubility score 0.594, disulfide risk True – use standard inclusion body conditions to maximise yield; refolding will follow.*

1. Inoculate 1 L LB (50 µg/mL kanamycin) in a 2.5 L baffled flask with 10 mL pre‑culture.
2. Incubate at 37 °C, 200 rpm until OD₆₀₀ reaches 0.6–0.8.
3. Add IPTG to a final concentration of **1.0 mM** (stock: 1 M in water, sterile‑filtered).
4. Continue incubation at **37 °C** for **4 h** with shaking (200 rpm).
5. For borderline solubility, an alternative low‑temperature induction (25 °C, 6 h, 0.2 mM IPTG) may be used, but the present protocol favours high inclusion body yield.

## 5. Cell Harvest and Lysis
- **Denaturing lysis is required** (warning: 6 cysteine residues; cytoplasmic disulfide formation is inefficient) (PMID:38479588).

1. Harvest cells by centrifugation at 6,000 × g, 4 °C, 15 min.
2. Wash pellet with 100 mL ice‑cold PBS (pH 7.4). Centrifuge again.
3. Resuspend pellet in **Denaturing Lysis Buffer** (6 M guanidine‑HCl, 20 mM Tris‑HCl pH 8.0, 500 mM NaCl, 20 mM imidazole) at 10 mL per gram wet pellet.
4. Stir at room temperature for 1 h (or until homogeneous). Optionally sonicate (6 × 10 s pulses, 30 % amplitude) to reduce viscosity.
5. Clarify lysate by centrifugation at 20,000 × g, 4 °C, 30 min. Keep supernatant (contains solubilised proinsulin).

## 6. Affinity Purification (IMAC) under Denaturing Conditions
1. Equilibrate 5 mL Ni‑NTA agarose (Qiagen) with **Denaturing Lysis Buffer** (10 column volumes, CV).
2. Load clarified lysate onto column at 1 mL/min.
3. Wash with 20 CV of **Wash Buffer** (6 M guanidine‑HCl, 20 mM Tris‑HCl pH 8.0, 500 mM NaCl, 50 mM imidazole).
4. Elute with **Elution Buffer** (6 M guanidine‑HCl, 20 mM Tris‑HCl pH 8.0, 500 mM NaCl, 500 mM imidazole). Collect 1 mL fractions.
5. Pool fractions with A₂₈₀ > 0.1. Analyse by SDS‑PAGE (expected ~13 kDa including tag).

## 7. Tag Removal
The 6×His‑tag is cleaved by Enterokinase after refolding (EK cleaves the sequence DDDDK†). Cleavage before refolding may expose hydrophobic surfaces and reduce yields. Therefore tag removal is performed **after refolding** (Section 8).

After refolding (Section 8), proceed:
1. Adjust refolded protein buffer to 20 mM Tris‑HCl pH 8.0, 50 mM NaCl, 2 mM CaCl₂. Do not include EDTA.
2. Add Enterokinase (EK) at a 1:100 (w/w) enzyme‑to‑substrate ratio.
3. Incubate at 16 °C for 16 h (or room temperature for 4 h). Confirm cleavage by SDS‑PAGE (shift from ~13 kDa to ~10 kDa).
4. Remove EK and uncleaved protein by passing through a 1 mL Ni‑NTA column (equilibrated in 20 mM Tris pH 8.0, 50 mM NaCl). Flow‑through contains tag‑free proinsulin.
5. Concentrate to ~1 mg/mL using a 3 kDa MWCO centrifugal filter.

## 8. Inclusion Body Refolding
- **Required because disulfide_risk = True** (6 cysteines must form three correct disulfide bonds).

### 8.1. On‑column refolding (preferred) 
1. After IMAC elution (Step 6), dilute the pooled eluate to 0.5–1 mg/mL protein in **Denaturing Refolding Buffer** (6 M urea, 20 mM Tris‑HCl pH 8.5, 500 mM NaCl, 1 mM reduced glutathione (GSH), 0.1 mM oxidised glutathione (GSSG), 0.5 M arginine).
2. Load onto a fresh Ni‑NTA column (equilibrated with the same buffer).
3. Wash with 5 CV of **Refolding Buffer** containing decreasing urea gradient (6 M → 2 M urea, stepwise 2 M per wash, 3 CV each).
4. Wash with 5 CV of **Final Refolding Buffer** (2 M urea, 20 mM Tris‑HCl pH 8.5, 500 mM NaCl, 1 mM GSH, 0.1 mM GSSG, 0.5 M arginine).
5. Incubate on‑column at 4 °C for 12–16 h.
6. Elute with **Elution Buffer** (same as step 4, but urea‑free: 20 mM Tris‑HCl pH 8.0, 500 mM NaCl, 500 mM imidazole). Collect fractions.

### 8.2. Alternative: dialysis refolding
1. Dialyse IMAC eluate (diluted to 0.2–0.5 mg/mL) against **Refolding Buffer** (without urea) containing 1 mM GSH, 0.1 mM GSSG, 0.5 M arginine, 20 mM Tris‑HCl pH 8.5, 50 mM NaCl. Use stepwise dialysis: 4 M urea → 2 M urea → 0 M urea, each step 4 h at 4 °C.
2. Oxidative refolding occurs during the final dialysis (18 h at 4 °C).

### 8.3. Post‑refolding cleanup
- Centrifuge refolded sample at 20,000 × g, 4 °C, 20 min to remove aggregates.
- Proceed to **Tag Removal** (Section 7) or directly to quality control.

## 9. Quality Control
1. **SDS‑PAGE** (15 % gel) under reducing and non‑reducing conditions. Expected mass: tag‑free proinsulin ~10 kDa. Non‑reducing band should be monomeric (no disulfide‑linked oligomers).
2. **Western blot** using anti‑insulin antibody (e.g., mouse monoclonal, 1:1000) to confirm identity.
3. **Analytical RP‑HPLC** (C18 column, 0.1 % TFA, acetonitrile gradient) – proinsulin elutes as a sharp peak; presence of isomers indicates misfolding.
4. **Mass spectrometry** (ESI‑MS) – expected mass 9393 Da (proinsulin without tag). Deviation > 1 Da suggests incomplete disulfide formation.
5. **Disulfide mapping** (peptide mapping with trypsin, LC‑MS/MS) to verify three disulfide bonds (Cys⁶–Cys¹¹, Cys⁷–Cys⁷?, actually human insulin: A6–A11, A7–B7, A20–B19; for proinsulin check appropriate residues).

## 10. Expected Yields and Troubleshooting Tips
- **Expected yield**: 20–50 mg of purified proinsulin per litre of culture after refolding and cleavage. Recovery from inclusion bodies is typically 10–30 % of total expressed protein.
- **Low refolding yield**: Optimise the GSH:GSSG ratio (commonly 10:1 to 3:1). Increase arginine to 0.8 M. Lower protein concentration during refolding (< 0.5 mg/mL).
- **Aggregation after cleavage**: Add 0.1 % Tween‑20 or increase NaCl to 300 mM. Perform cleavage at 4 °C.
- **Proteolytic degradation**: Include 1 mM PMSF during lysis; use protease‑inhibitor cocktail. Ensure all buffers are ice‑cold.
- **Incomplete tag cleavage**: Increase EK amount or incubation time. Check Ca²⁺ concentration (2 mM required).
- **Expression**: If no induction observed, sequence the plasmid; ensure correct reading frame. Can co‑express GroEL/GroES (pGro7) to reduce aggregation (optional, not mandatory for inclusion body route).
- **Note**: To obtain mature insulin, the C‑peptide must be enzymatically removed from proinsulin (e.g., trypsin and carboxypeptidase B). This step is **not** included in the present protocol.