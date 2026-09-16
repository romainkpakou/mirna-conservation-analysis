#!/usr/bin/env python3
"""
Croise la conservation observée par alignement de séquence (script 02) avec
l'orthologie annotée dans Ensembl Compara -- confirme que la conservation
mesurée reflète bien des orthologues réels, pas une similarité de séquence
fortuite.

Limité aux miARN humains à locus génomique unique et non ambigu (MIR208A,
MIR208B, MIR499A) : miR-1, miR-133a et let-7a existent chez l'humain sous
plusieurs copies génomiques paralogues (ex: MIR1-1/MIR1-2), ce qui rend le
choix d'un gène Ensembl "canonique" arbitraire pour ces trois familles --
non traité ici plutôt que de trancher arbitrairement (voir PLAN.md).
"""
import csv
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

ENSEMBL_BASE = "https://rest.ensembl.org"

# Gènes humains à locus unique, vérifiés par coordonnées GRCh38 (voir PLAN.md)
FAMILY_TO_GENE = {
    "miR-208a-3p": "ENSG00000199157",  # MIR208A, intronique dans MYH6
    "miR-208b-3p": "ENSG00000215991",  # MIR208B, intronique dans MYH7
    "miR-499a-5p": "ENSG00000207635",  # MIR499A, intronique dans MYH7B
}

TARGET_SPECIES = ["mouse", "rattus_norvegicus", "chicken", "zebrafish"]


def fetch_homologies(gene_id: str) -> list[dict]:
    all_homologies = []
    for target in TARGET_SPECIES:
        url = f"{ENSEMBL_BASE}/homology/id/human/{gene_id}?content-type=application/json;target_species={target}"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read())
            for entry in data.get("data", []):
                for h in entry.get("homologies", []):
                    all_homologies.append({
                        "target_species": h["target"]["species"],
                        "homology_type": h["type"],
                        "pct_identity_gene_level": h["target"]["perc_id"],
                        "taxonomy_level": h.get("taxonomy_level", "NA"),
                    })
        except Exception as e:
            print(f"  ! échec requête Compara pour {target}: {e}")
    return all_homologies


def main():
    rows = []
    for family, gene_id in FAMILY_TO_GENE.items():
        print(f"{family} ({gene_id}) : interrogation Ensembl Compara...")
        homologies = fetch_homologies(gene_id)
        for h in homologies:
            rows.append({"famille": family, "ensembl_gene_id": gene_id, **h})
        print(f"  {len(homologies)} orthologue(s) trouvé(s) : "
              + ", ".join(f"{h['target_species']} ({h['homology_type']})" for h in homologies))

    out_path = RESULTS / "compara_orthology_crosscheck.tsv"
    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=[
            "famille", "ensembl_gene_id", "target_species", "homology_type",
            "pct_identity_gene_level", "taxonomy_level"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nTerminé -> {out_path}")


if __name__ == "__main__":
    main()
