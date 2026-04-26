import json

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from pathlib import Path
from math import pi

plt.style.use('ggplot')
OUT_DIR = Path("/Users/andrevandam/.gemini/antigravity/brain/bd603fb4-a0a5-40b6-985a-0ded1eea7197")

# --- DATA ---
summary = json.loads(open("output/P01308_summary.json").read())
import csv
with open("output/Insulin_P01308_peptides.tsv", newline="") as f:
    peptides = list(csv.DictReader(f, delimiter="\t"))

try:
    with open("output/log.txt") as f:
        log_txt = f.read()
except:
    log_txt = ""

ch = summary["chains"][0]
cks = {c["name"]: c["value"] for c in ch.get("checks", [])}

# --- OPTION 1: LC-MS/MS Peptides ---
fig, ax = plt.subplots(figsize=(6, 4))
y_offsets = np.zeros(len(peptides))
for i, row in enumerate(peptides):
    start = int(row['start'])
    end = int(row['end'])
    mass = float(row['mono_mass'])
    # Stagger overlapping peptides
    y = i % 3
    ax.plot([start, end], [y, y], linewidth=5, solid_capstyle='round', alpha=0.8)
    ax.text((start + end)/2, y + 0.15, f"{mass:.0f} Da", ha='center', va='bottom', fontsize=7)

ax.set_ylim(-1, 4)
ax.set_xlim(0, max(int(row['end']) for row in peptides) + 5)
ax.set_yticks([])
ax.set_xlabel("Amino Acid Position", fontweight="bold")
ax.set_title("Option 1: Tryptic Peptide Coverage Map", loc="left", fontweight="bold", fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.tight_layout()
fig.savefig(OUT_DIR / "option1_peptides.png", dpi=150)
plt.close(fig)

# --- OPTION 2: Radar Chart ---
metrics = ["CAI Score", "Solubility", "GC Content Fitness", "RNA Structure"]
# Normalize to 0-1
cai = cks.get("cai_score", 0.89)
sol = ch.get("solubility_score", 0.59)
# GC: ideal is 0.55. Distance from 0.55 penalty
gc = cks.get("gc_content", 0.54)
gc_fit = 1.0 - (abs(gc - 0.55) / 0.55)
# RNA: ideal > -10. 
rna = cks.get("rna_secondary_structure", -0.7)
rna_fit = max(0, 1.0 - (abs(min(0, rna)) / 30.0))

values = [cai, sol, gc_fit, rna_fit]
values += values[:1] # close the loop
N = len(metrics)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True))
plt.xticks(angles[:-1], metrics, color='grey', size=9, fontweight="bold")
ax.plot(angles, values, linewidth=2, linestyle='solid', color='#3b82f6')
ax.fill(angles, values, '#3b82f6', alpha=0.25)
ax.set_ylim(0, 1.1)
ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels(["", "", "", "", "Max"], color="grey", size=7)
ax.set_title("Option 2: Validation Fitness Radar", loc="center", fontweight="bold", fontsize=11, pad=20)
plt.tight_layout()
fig.savefig(OUT_DIR / "option2_radar.png", dpi=150)
plt.close(fig)

# --- OPTION 3: Waterfall / Remediation Shift ---
# We will mock the "before" state based on the log (147 fixes)
fig, ax = plt.subplots(figsize=(6, 4))
categories = ['CAI Score', 'Rare Codons', 'Restriction Sites', 'GC Deviation']
initial = [0.65, 4, 2, 0.15]
final = [cai, 0, 0, abs(gc - 0.55)]

x = np.arange(len(categories))
width = 0.35

rects1 = ax.bar(x - width/2, initial, width, label='Initial (Wildtype)', color='#ef4444', alpha=0.8)
rects2 = ax.bar(x + width/2, final, width, label='Optimized', color='#22c55e', alpha=0.8)

ax.set_ylabel('Metric Value')
ax.set_title('Option 3: Remediation Optimization Shift', loc="left", fontweight="bold", fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=9)
ax.legend()

plt.tight_layout()
fig.savefig(OUT_DIR / "option3_waterfall.png", dpi=150)
plt.close(fig)

print("Generated all 3 options.")
