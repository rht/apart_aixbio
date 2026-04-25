"""Visualize a CIF structure file as a static 3D CA-backbone plot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser


def extract_ca_coords(cif_path: str) -> np.ndarray:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(Path(cif_path).stem, cif_path)
    coords = [
        residue["CA"].get_vector().get_array()
        for model in structure
        for chain in model
        for residue in chain
        if residue.id[0] == " " and "CA" in residue
    ]
    if not coords:
        print(f"No CA atoms found in {cif_path}", file=sys.stderr)
        sys.exit(1)
    return np.array(coords)


def plot_backbone(coords: np.ndarray, title: str, out_path: str | None) -> None:
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(*coords.T, linewidth=1.5, color="steelblue")
    ax.scatter(*coords[0], color="green", s=40, label="N-term")
    ax.scatter(*coords[-1], color="red", s=40, label="C-term")
    ax.set_title(title)
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.legend()
    if out_path:
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {out_path}")
    else:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", help="Path to a .cif structure file")
    parser.add_argument("-o", "--output", help="Save plot to file instead of displaying")
    parser.add_argument("-t", "--title", help="Plot title (default: filename stem)")
    args = parser.parse_args()

    coords = extract_ca_coords(args.cif)
    title = args.title or Path(args.cif).stem
    plot_backbone(coords, title, args.output)


if __name__ == "__main__":
    main()
