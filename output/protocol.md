# Protocol: Insulin (Proinsulin) Expression and Purification in *E. coli*

## Literature Context
Human insulin was the first recombinant therapeutic protein produced in *E. coli* (PMID:85300) by expressing separate A and B chains fused to β‑galactosidase, followed by chemical cleavage and oxidative refolding. Modern strategies often express proinsulin as a single chain (e.g., using strain *E. coli* 20; PMID:30735706). *E. coli* remains the most cost‑effective host for insulin production, although the reducing cytoplasm prevents proper disulfide formation, necessitating denaturing purification and refolding (PMID:38479588, PMID:25270715). For the proinsulin sequence (86 aa, 6 cysteines) with a 6×His tag, the predicted solubility (score 0.594) is high, but the presence of three disulfide bonds mandates inclusion body processing and oxidative refolding.

## 1. Strains, Plasmid, and Safety Notes
- **Host strain:** *E. coli* BL21‑CodonPlus‑RIL (DE3) or BL21(DE3)Star (for high‑level T7‑driven expression).
- **Plasmid:** pET‑28a(+) encoding **6×His‑Enterokinase‑proinsulin** (Insulin P01308, aa 25–110) cloned via *BamHI/XhoI*.
- **Safety:** BSL‑1. All waste containing imidazole or β‑mercaptoethanol must be collected and disposed according to institutional guidelines.
- **Key reagents:** Kanamycin (50 µg/mL), IPTG, urea (ultrapure), glutathione (GSH/GSSG), enterokinase (Ek, e.g., GenScript E‑165‑1).

## 2. Transformation
1. Thaw 50 µL electrocompetent BL21(DE3) cells on ice.
2. Add 1 µL (10–100 ng) plasmid DNA, transfer to a 0.1 cm cuvette.
3. Electroporate at 1.8 kV, 25 µF, 200 Ω.
4. Immediately add 1 mL SOC medium, incubate 1 h at 37 °C.
5. Spread 50–200 µL on LB‑kanamycin (50 µg/mL) agar; incubate overnight at 37 °C.

## 3. Pre‑culture
- Pick a single colony into 5 mL LB‑kan (50 µg/mL).
- Grow at 37 °C, 220 rpm, until OD₆₀₀ ≈ 0.6–0.8 (≈4–5 h).
- Use 1 % (v/v) to inoculate expression culture.

## 4. Expression Culture and Induction
- **Medium:** 1 L LB‑kan (50 µg/mL) in a 2.5 L baffled flask.
- **Growth:** 37 °C, 220 rpm until OD₆₀₀ = 0.6–0.8.
- **Induction:** Add IPTG to **0.5 mM** (optimal for solubility score >0.45; PMID:38479588).
- **Conditions:** 37 °C for **4 h** (or 25 °C overnight if lower inclusion body yield is tolerable).
- **Note:** Predicted solubility is high, but due to 6 cysteines the protein will accumulate in inclusion bodies. Induction at 37 °C maximises inclusion body yield without harming cell growth.

## 5. Cell Harvest and Lysis
1. Harvest by centrifugation at 6 000 × g, 4 °C, 15 min.
2. Discard supernatant, resuspend pellet in **30 mL Lysis Buffer** (50 mM Tris‑HCl pH 8.0, 150 mM NaCl, 1 mM EDTA, 1 mM PMSF, 0.1 % Triton X‑100).
3. Lyse by sonication (6 × 30 s pulses, 40 % amplitude, on ice).
4. Centrifuge lysate at 20 000 × g, 4 °C, 30 min. Discard supernatant (contains soluble host proteins).
5. **Wash inclusion bodies:** resuspend pellet in 30 mL **Wash Buffer 1** (50 mM Tris‑HCl pH 8.0, 1 M urea, 2 % Triton X‑100), sonicate briefly, centrifuge (20 000 × g, 15 min). Repeat once.
6. Wash twice with **Wash Buffer 2** (50 mM Tris‑HCl pH 8.0, 1 M urea, no detergent) to remove residual detergent.

## 6. Affinity Purification (IMAC) under Denaturing Conditions
- **Denaturing buffer:** 50 mM Tris‑HCl pH 8.0, 8 M urea, 300 mM NaCl, 20 mM imidazole, 1 mM DTT (added fresh).
- **Solubilization:** Resuspend washed inclusion bodies in 20 mL denaturing buffer. Stir at RT for 1 h; centrifuge at 20 000 × g, 30 min (RT). Keep supernatant.
- **Column:** 5 mL Ni‑NTA agarose (GE Healthcare) equilibrated with denaturing buffer.
- **Load:** Apply clarified solubilised protein at 1 mL/min (RT).
- **Wash:** 10 column volumes (CV) of denaturing buffer + 40 mM imidazole.
- **Elution:** 5 CV denaturing buffer + 300 mM imidazole. Collect 1 mL fractions. Pool those containing proinsulin (A₂₈₀ monitoring).

## 7. Tag Removal (Enterokinase Cleavage)
- **Buffer exchange:** Dialyse pooled IMAC eluate against **EcoBuffer** (50 mM Tris‑HCl pH 8.0, 150 mM NaCl, 2 mM CaCl₂) containing 2 M urea (to maintain solubility) at 4 °C, overnight, using 10 kDa MWCO membrane.
- **Digestion:** Add Ek at 1 U per 100 µg target protein, incubate at 25 °C for **18 h** (or follow manufacturer’s recommendation).
- **Separation:** Pass through a fresh Ni‑NTA column (0.5 mL resin) equilibrated with EcoBuffer + 2 M urea. Collect flow‑through (untagged proinsulin). Wash with 2 CV of same buffer. Pool flow‑through and wash – this removes the 6×His‑Ek fragment.

## 8. Inclusion Body Refolding (Oxidative Refolding of Proinsulin)
> **Required because disulfide_risk = True.**  
> The protein contains 6 cysteines that must form three native disulfides (A6–A11, A7–B7, A20–B19). Refolding is performed on the fully reduced, denatured proinsulin (obtained after tag removal).

1. **Adjust protein concentration** to **0.2–0.5 mg/mL** (A₂₈₀, ε = 1.2 mL mg⁻¹ cm⁻¹ predicted) in denaturing buffer (8 M urea, 50 mM Tris‑HCl pH 8.0, 1 mM DTT).  
2. **Refolding buffer:** 50 mM Tris‑HCl pH 8.0, 1 M L‑arginine, 2 mM EDTA, **1 mM reduced glutathione (GSH)**, **0.1 mM oxidised glutathione (GSSG)** (redox ratio 10:1, optimised for insulin; PMID:25270715).  
3. **Rapid dilution:** Slowly add denatured protein (with gentle stirring) into 10 volumes of refolding buffer at 4 °C. Final urea concentration ≤ 0.8 M.  
4. **Incubation:** Stir at 4 °C for **24–48 h** in the dark. Monitor turbidity; if precipitation occurs, dilute further or lower temperature to 10 °C.  
5. **Dialysis (optional):** Dialyse refolded protein against 50 mM Tris‑HCl pH 8.0, 100 mM NaCl (no argon) at 4 °C to remove arginine and glutathione.  
6. **Concentration:** Use centrifugal concentrator (3 kDa MWCO) to achieve 0.5–1.0 mg/mL.

## 9. Quality Control
- **SDS‑PAGE** (15 % gel, reducing vs. non‑reducing): Reduced proinsulin runs as ~12 kDa; correctly formed disulfides cause a mobility shift (higher apparent mass under non‑reducing).  
- **Analytical RP‑HPLC** (C18 column, 0.1 % TFA in water/acetonitrile): Retention time relative to proinsulin standard; purity ≥ 95 %.  
- **Mass spectrometry** (ESI‑MS): Expected mass of untagged proinsulin = 9 396 Da (reduced), 9 390 Da (oxidised, all disulfides formed).  
- **Insulin activity** (optional): After C‑peptide removal (trypsin + carboxypeptidase B; not covered in this protocol), perform insulin receptor binding assay or glucose uptake assay.  
- **Endotoxin**: < 1 EU/mg if intended for cell‑based assays.

## 10. Expected Yields and Troubleshooting Tips
- **Expected yield:** 30–60 mg of purified proinsulin per litre of culture (after tag removal and refolding).  
- **Low inclusion body yield:** Ensure IPTG is added at OD₆₀₀ 0.6–0.8; use fresh transformants.  
- **Refolding precipitation:** Reduce protein concentration to 0.1 mg/mL; increase GSH/GSSG ratio to 5:1 or add 0.5 M L‑arginine; lower temperature to 10 °C.  
- **Incomplete disulfide formation:** Supplement redox pair with 0.2 mM cystine; extend incubation to 72 h; degas buffers with N₂.  
- **Tag removal failure:** Confirm Ek activity; ensure Ca²⁺ ≥ 1 mM; include 0.1 % Tween‑20 in digestion buffer.  
- **Final insulin production:** The protocol yields refolded proinsulin. To obtain mature insulin, digest with trypsin (1:100, w/w) and carboxypeptidase B (1:500, w/w) at 25 °C for 30 min (stop with PMSF). Purify by preparative RP‑HPLC.

**Literature References**  
- PMID:85300 – First expression of insulin A and B chains in *E. coli*.  
- PMID:33340315 – Historical perspective on Humulin development.  
- PMID:30735706 – Proinsulin expression in *E. coli* strain 20.  
- PMID:38479588 – General challenges in *E. coli* recombinant protein production.  
- PMID:25270715 – Cell factories and refolding strategies for insulin.